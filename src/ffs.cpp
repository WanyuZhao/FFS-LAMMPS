#include<unistd.h>
#include<stdarg.h>
#include<stdio.h>
#include<string.h>
#include<mpi.h>
#include<cstdlib>
#include<ctime>
#include<map>
#include<queue>
#include<string>
#include"lammps.h"
#include"library.h"
#include"input.h"
#include"update.h"
#include"ffs.h"
using namespace LAMMPS_NS;
#define DEBUG printf("------ rank %d (%d of universe %d) ------ line %d ------\n", world->rank, local->rank, local->id, __LINE__);
//the class that contains the information of process
struct MpiInfo {
    const int LEADER;
    //id is equal to myShooting
    const int id;
    MPI_Comm comm;
    int size,rank;
    bool isLeader;
    MpiInfo(MPI_Comm comm,int id=0):comm(comm),LEADER(0),id(id) {
        MPI_Comm_size(comm,&size);
        MPI_Comm_rank(comm,&rank);
        //if the rank is 0, the isLeader is true
        isLeader=LEADER==rank;
    };
};
//local is the communicator for lammps, each process in local runs a lammps instance
//world is the communicator for ffs, each process in world runs a ffs
const MpiInfo *world,*local;
//The class is used for parallel execution
class FfsBranch {
public:
    FfsBranch() {
      if (!FfsBranch::commInited) {
        //it's equal to the shooting number
        size = world->size / local->size;
        FfsBranch::initComm();
        FfsBranch::commInited=true;
      }
    }
    ~FfsBranch() {
    }
private:
    static bool commInited;
    static void initComm() {
        MPI_Comm_dup(local->comm, &FfsBranch::commLocal);
        //split the leader process from local, and assign it to the commLeader
        MPI_Comm_split(world->comm, local->isLeader ? 0 : MPI_UNDEFINED, world->rank, &FfsBranch::commLeader);
    };
protected:
    static int size;
    static MPI_Comm commLeader,commLocal;
    static const int TAG_COUNTDOWN_DONE=1;
    static const int TAG_COUNTDOWN_TERMINATE=2;
    static const int TAG_FILEWRITER_LINE=3;
    static const int TAG_FILEREADER=5;
    static const int TAG_STATS_FLUSH=6;
};
bool FfsBranch::commInited=false;
//construct 2 communicator, commLeader contains the leader process, 
MPI_Comm FfsBranch::commLeader,FfsBranch::commLocal;
//it's equal to the cpuEach
int FfsBranch::size=0;
struct FfsFileReader: public FfsBranch {
    std::map<std::string,std::string> dict;
    //Initialization function, with input parameter "filename" in char pointer type. It can construct a dictionary with "key" refers to the command, and "value" refers to the parameter
    FfsFileReader(const char *filename) {
        //only the world leader process read file
        if (world->isLeader) {
            FILE *f=fopen(filename,"r");
            static char cLine[1024];
            while (fgets(cLine,sizeof(cLine),f)) {
                std::string sLine(cLine);
                int l=sLine.find_first_of('#');
                //check if "#" exist. If true, substract the substring for the beginning to the "#" such as "command    parameter   #comment"
                if (l!=std::string::npos) {
                    sLine=sLine.substr(0,l);
                }
                const std::string SPACE=" \t\n\v\f\r";
                int begin1=sLine.find_first_not_of(SPACE);
                //delete the white space characters
                if (begin1==std::string::npos) {
                    continue;
                }
                int end1=sLine.find_first_of(SPACE,begin1);
                //find the end of command, substract the command, like "command parameter"
                if (end1==std::string::npos) {
                    end1=sLine.length();
                }
                std::string key=sLine.substr(begin1,end1-begin1);
                dict[key]="";
                //substract the parameter
                int begin2=sLine.find_first_not_of(SPACE,end1);
                if (begin2==std::string::npos) {
                    continue;
                }
                int end2=sLine.find_last_not_of(SPACE)+1;
                std::string value=sLine.substr(begin2,end2-begin2);
                dict[key]=value;
            }
            fclose(f);
        }
    };

    //Parses a command string "parameter" into a int value
    int getInt(const std::string &name) const {
        int y=0;
        const std::string v=getString(name);
        sscanf(v.c_str(),"%d",&y);
        return y;
    }

    //Parses a command string "<param1, param2, ...>" into a vector of parameters
    std::vector<int> getVector(const std::string &name) const {
        std::vector<int> v;
        const std::string s = getString(name);
        const char *p = s.c_str();
        int x;
        int n;
        while (sscanf(p, "%d%n", &x, &n) > 0) {
            v.push_back(x);
            p += n;
        }
        return v;
    }
    //The function is used for get the command parameter and, it can be executed in parallel environment. Input value is the string of command
    const std::string getString(const std::string &name) const {
        char *p;
        int l;
        if (world->isLeader) {
            //get the value of the key
            std::map<std::string,std::string>::const_iterator i=dict.find(name);
            if (i==dict.end()) {
              fprintf(stderr, "Missing parameter \"%s\" in ffs input\n", name.c_str());
              p = NULL;
              l = 0;
            }
            else {
                const std::string s = i->second;
                l = s.length();
                p = new char[l + 1];
                s.copy(p, l);
            }
        }
        //here for parallel excution
        MPI_Bcast(&l, 1, MPI_INT, 0, world->comm);
        if (!world -> isLeader) {
            p = new char[l + 1];
        }
        MPI_Bcast(p, l + 1, MPI_CHAR, 0, world->comm);
        p[l] = '\0';
        const std::string result = std::string(p);
        delete[] p;
        return result;
    }
};
const FfsFileReader *ffsParams;
class FfsFileWriter: public FfsBranch {
protected:
    FfsFileWriter(const char *filename) {
        if (world->isLeader) {
            f=fopen(filename,"w");
        }
        nFlush=0;
    }
    ~FfsFileWriter() {
        check();
        if (world->isLeader) {
            fclose(f);
        }
    }

    //receive a series of string, and write it into file
    void writeln(const char *format,...) {
        if (!local->isLeader) {
            return ;
        }
        static char buffer[MAX_LENGTH+1];
        //variable list
        va_list args;
        //format the input and write it into buffer
        va_start(args,format);
        vsprintf(buffer,format,args);
        va_end(args);
        //if it's main process, write it into file, if not, send it to the main process
        putstr(buffer);
    }

    //check if there exists information from the non-main process, if true, receive it and write it into file
    void check() {
        if (!world->isLeader) {
            return ;
        }
        static char buffer[MAX_LENGTH];
        //continuosly receive the broadcast from the process which is not the main process, and then write it into file
        while (1) {
            int flag;
            MPI_Status status;
            MPI_Iprobe(MPI_ANY_SOURCE, FfsBranch::TAG_FILEWRITER_LINE, FfsBranch::commLeader, &flag, &status);
            if (flag) {
                int l;
                MPI_Get_count(&status,MPI_CHAR,&l);
                printf("[date=%d] world leader will receive TAG_FILEWRITER_LINE from %d\n", std::time(0), status.MPI_SOURCE);
                //receive the information and store it in the buffer
                MPI_Recv(buffer, l, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, FfsBranch::commLeader, &status);
                printf("[date=%d] world leader did receive TAG_FILEWRITER_LINE from %d\n", std::time(0), status.MPI_SOURCE);
                putstr0(status.MPI_SOURCE,buffer);
            }
            else {
                break;
            }
        }
    }

    //tell that there's n more lines in the buffer
    void pushGroup(int n) {
        nFlush+=n;
    }
private:
    FILE *f;
    //if the process is not main process, the function will send the char s to the main process
    void putstr(char *s) {
        s[MAX_LENGTH-1]='\0';
        if (world->isLeader) {
            putstr0(0,s);
        }
        else {
            //the length does not include the terminator
            int l=strlen(s);
            printf("[date=%d] universe %d will send TAG_FILEWRITER_LINE\n", std::time(0), local->id);
            //send the value
            MPI_Send(s, l + 1, MPI_CHAR, 0, FfsBranch::TAG_FILEWRITER_LINE, FfsBranch::commLeader);
            printf("[date=%d] universe %d did send TAG_FILEWRITER_LINE\n", std::time(0), local->id);
        }
    }
    int nFlush;

    //write the sender and char s into the buffer, and then write the buffer into file, nflush stands for the number of lines need to be written
    void putstr0(int sender,const char *s) {
        static char buffer[MAX_LENGTH];
        sprintf(buffer,"%4d    %s",sender,s);
        fprintf(f,"%s\n",buffer);
        printf("%s\n",buffer);
        if (nFlush<=0) {
            fflush(f);
        }
        else {
            nFlush--;
        }
    }
    static const int MAX_LENGTH=100;
};
class FfsLambdaLogger: public FfsFileWriter {
public:
    FfsLambdaLogger():FfsFileWriter("lambda.txt") {
    }
    void check() {
        FfsFileWriter::check();
    }

    //write the xyzInit as heading, xyzFinal as ending, and the middle each line is a element of the vector
    void writeln(const char *xyzInit, const char *xyzFinal,const std::vector<int> &v) {
        int n=v.size();
        //1 for heading, 1 for ending, and n for n elements in vector
        FfsFileWriter::pushGroup(1+n+1);
        FfsFileWriter::writeln("BEGIN %s",xyzInit);
        int i;
        for (i=0;i<n;i++) {
            FfsFileWriter::writeln("%d",v[i]);
        }
        FfsFileWriter::writeln("END %s",xyzFinal);
    }
};
class FfsTrajectoryWriter: public FfsFileWriter {
public:
    FfsTrajectoryWriter():FfsFileWriter("trajectory.out.txt") {
    }
    void check() {
        FfsFileWriter::check();
    }
    //write the parameters into the "trajectory.out.txt"
    void writeln(const char *xyzInit,int lambdaInit,int velocitySeed,int64_t timestep,const char *xyzFinal,int lambdaFinal) {
        if (xyzInit==0) {
            FfsFileWriter::writeln("___ (__________)  >==%010d %20lld==>  %3d (xyz.%s)",velocitySeed,timestep,lambdaFinal,xyzFinal);
        }
        else {
            FfsFileWriter::writeln("%3d (xyz.%s)  >==%010d %20lld==>  %3d (xyz.%s)",lambdaInit,xyzInit,velocitySeed,timestep,lambdaFinal,xyzFinal);
        }
    }
};
class FfsTrajectoryReader: FfsBranch {
public:
    FfsTrajectoryReader() {
        int n;
        int *p;
        if (world->isLeader) {
            std::vector< std::vector<int> > v;
            v.resize(FfsBranch::size);
            FILE *f=fopen("trajectory.in.txt","r");
            //attention: here only the world part do the while, v is a vector with structure as v[brunch] = {layer, count, lambda}
            while (1) {
                int lambda,layer,branch,count;
                //ignore the string and char, read the int number and store in the 4 variable
                int ret=fscanf(f,"%*s%*s%*s%*s%*s%d%*[^.]%*c%d%*c%*c%d%*c%d%*[^\n]%*c",&lambda,&layer,&branch,&count);
                if (ret==EOF) {
                    break;
                }
                //it's a two-dimensional vector
                v[branch].push_back(layer);
                v[branch].push_back(count);
                v[branch].push_back(lambda);
            }
            fclose(f);
            int i,j;
            if (1) {
                //copy the first vector in v to the pointer "p" address
                const std::vector<int> vv=v[0];
                n=vv.size();
                p=new int[n];
                std::copy(vv.begin(),vv.end(),p);
            }
            for (int i=1;i<(int)v.size();i++) {
                //copy the information and send it to the local leader process
                const std::vector<int> &vv=v[i];
                int nn=vv.size();
                int *pp=new int[nn];
                for (j=0;j<nn;j++) {
                    pp[j]=vv[j];
                }
                printf("[date=%d] world leader will send TAG_FILEREADER to %d\n", std::time(0), i);
                MPI_Send(pp, nn, MPI_INT, i, FfsBranch::TAG_FILEREADER, FfsBranch::commLeader);
                printf("[date=%d] world leader did send TAG_FILEREADER to %d\n", std::time(0), i);
                delete[] pp;
            }
        }
        else if (local->isLeader) {
            MPI_Status status;
            //check if there is information in the communicator, and count it
            MPI_Probe(0, FfsBranch::TAG_FILEREADER, FfsBranch::commLeader, &status);
            MPI_Get_count(&status,MPI_INT,&n);
            p=new int[n];
            printf("[date=%d] universe %d will receive TAG_FILEREADER\n", std::time(0), local->id);
            //receive and store it in the pointer p
            MPI_Recv(p, n, MPI_INT, 0, FfsBranch::TAG_FILEREADER, FfsBranch::commLeader, &status);
            printf("[date=%d] universe %d did receive TAG_FILEREADER\n", std::time(0), local->id);
        }
        //receive and store value
        if (local->isLeader) {
            int i;
            for (i=0;i+2<n;i+=3) {
                int layer=p[i];
                int count=p[i+1];
                int lambda=p[i+2];
                if (layer+1>lambdaLocal.size()) {
                    lambdaLocal.resize(layer+1);
                }
                std::vector<int> &v=lambdaLocal[layer];
                if (count+1>v.size()) {
                    v.resize(count+1);
                }
                v[count]=lambda;
            }
            delete[] p;
        }
        //make sure all the process has been finished
        MPI_Barrier(FfsBranch::commLocal);
    }

    //get the vector of lambda according to the number of layer
    const std::vector<int> &get(int layer) const {
        if (layer<lambdaLocal.size()) {
            return lambdaLocal[layer];
        }
        else {
            return emptyVector;
        }
    }

    //count the total length of the variable, layer, in each process 
    int countPrecalculated(int layer) const {
        int x;
        if (local->isLeader) {
            int t=get(layer).size();
            //count the total length of layer
            MPI_Allreduce(&t, &x, 1, MPI_INT, MPI_SUM, FfsBranch::commLeader);
        }
        //send to all process
        MPI_Bcast(&x, 1, MPI_INT, 0, FfsBranch::commLocal);
        return x;
    }
private:
    std::vector< std::vector<int> > lambdaLocal;
    std::vector<int> emptyVector;
};
class FfsCountdown: public FfsBranch {
public:
    //when n <= 0, the terminated label turns to true
    FfsCountdown(int n) {
        remains=n;
        terminated=n<=0;
    }
    ~FfsCountdown() {
    }

    //send the variable n to the communicator
    //the done function must be continued by next() function
    void done(int n=1) {
        //only consider the local leader process
        if (!local->isLeader) {
            return ;
        }

        //check if the object has already been terminated
        if (terminated) {
            return ;
        }
        int x=n;
        if (world->isLeader) {
            remains-=x;
        }
        else {
            printf("[date=%d] universe %d will send TAG_COUNTDOWN_DONE\n", std::time(0), local->id);
            //send the variable, n, to the world leader process, with tag FfsBranch::TAG_COUNTDOWN_DONE
            MPI_Send(&x, 1, MPI_INT, 0, FfsBranch::TAG_COUNTDOWN_DONE, FfsBranch::commLeader);
            printf("[date=%d] universe %d did send TAG_COUNTDOWN_DONE\n", std::time(0), local->id);
        }
    }

    //check if there's next
    bool next() {
        if (terminated) {
            return false;
        }
        int ret;
        //only the local leader process handle the data
        if (local->isLeader) {
            if (world->isLeader) {
                //continuosly receive data
                while (1) {
                    int flag;
                    //store the MPI information
                    MPI_Status status;
                    //check if there's message tagged with countdown
                    MPI_Iprobe(MPI_ANY_SOURCE, FfsBranch::TAG_COUNTDOWN_DONE, FfsBranch::commLeader, &flag, &status);
                    if (flag) {
                        int x;
                        printf("[date=%d] world leader will receive TAG_COUNTDOWN_DONE from %d\n", std::time(0), status.MPI_SOURCE);
                        //receive the data from other process
                        MPI_Recv(&x, 1, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, FfsBranch::commLeader, &status);
                        printf("[date=%d] world leader did receive TAG_COUNTDOWN_DONE from %d\n", std::time(0), status.MPI_SOURCE);
                        remains-=x;
                    }
                    else {
                        break;
                    }
                }
                //there is no remains
                if (remains<=0) {
                    terminated=true;
                    int i;
                    for (i = 1; i < FfsBranch::size; i += 1) {
                        printf("[date=%d] world leader will send TAG_COUNTDOWN_TERMINATE to %d\n", std::time(0), i);
                        //send a message to all local leader for terminate
                        MPI_Send(0, 0, MPI_INT, i, FfsBranch::TAG_COUNTDOWN_TERMINATE, FfsBranch::commLeader);
                        printf("[date=%d] world leader did send TAG_COUNTDOWN_TERMINATE to %d\n", std::time(0), i);
                    }
                    ret=0;
                }
                else {
                    ret=1;
                }
            }
            else {
                int flag;
                MPI_Status status;
                //check if there's message about terminate
                MPI_Iprobe(0, FfsBranch::TAG_COUNTDOWN_TERMINATE, FfsBranch::commLeader, &flag, &status);
                if (flag) {
                    terminated=true;
                    printf("[date=%d] universe %d will receive TAG_COUNTDOWN_TERMINATE\n", std::time(0), local->id);
                    MPI_Recv(0, 0, MPI_INT, 0, FfsBranch::TAG_COUNTDOWN_TERMINATE, FfsBranch::commLeader, &status);
                    printf("[date=%d] universe %d did receive TAG_COUNTDOWN_TERMINATE\n", std::time(0), local->id);
                    ret=0;
                }
                else {
                    ret=1;
                }
            }
        }
        //broadcast the value of ret
        //let Local leader share the data to the non-leader process
        MPI_Bcast(&ret, 1, MPI_INT, 0, FfsBranch::commLocal);
        if (ret==0) {
            terminated=true;
        }
        return ret;
    }
private:
    //label the status of the object
    bool terminated;

    //the number of remains
    int remains;
};
class FfsFileTree: public FfsBranch {
public:
    FfsFileTree(const FfsTrajectoryReader *ftr,int layer):layer(layer) {
        a = new int[FfsBranch::size];
        b = new int[FfsBranch::size];
        int i;
        for (i = 0; i < FfsBranch::size; i += 1) {
            a[i]=0;
            b[i]=0;
        }
        if (local->isLeader) {
            //get a vector with components of the value of lambda
            lambdaLocal=ftr->get(layer);
            a[local->id]=lambdaLocal.size();
        }
    }
    ~FfsFileTree() {
        delete[] b;
        delete[] a;
    }

    //add a new lambda into the end of lambdaLocal vector, and return the new name
    const std::string add(int lambda) {
        int x;
        if (local->isLeader) {
            x=a[local->id];
            a[local->id]++;
            lambdaLocal.push_back(lambda);
        }
        MPI_Bcast(&x, 1, MPI_INT, 0, FfsBranch::commLocal);
        return generateName(x);
    }

    //Synchronize the values of lambda across all processes
    void commit() {
        if (local->isLeader) {
            //get the global value of lambda size and store in b
            MPI_Allreduce(a, b, FfsBranch::size, MPI_INT, MPI_SUM, FfsBranch::commLeader);
            int mySize=lambdaLocal.size();
            int currentSize;
            //store the total size of lambdaLocal in each process
            int maxSize;
            MPI_Allreduce(&mySize, &maxSize, 1, MPI_INT, MPI_MAX, FfsBranch::commLeader);
            int *p=new int[maxSize];
            int i,j;
            for (i = 0; i < FfsBranch::size; i += 1) {
                if (i==local->id) {
                    currentSize=mySize;
                    std::copy(lambdaLocal.begin(),lambdaLocal.end(),p);
                }
                //there is a deadlock, share the data to all the shooting
                MPI_Bcast(&currentSize, 1, MPI_INT, i, FfsBranch::commLeader);
                MPI_Bcast(p, currentSize, MPI_INT, i, FfsBranch::commLeader);
                //copy the p to the lambdaGlobal
                for (j=0;j<currentSize;j++) {
                    //lambdaGlobal stores the lambda value in all shooting
                    lambdaGlobal.push_back(p[j]);
                }
            }
            delete[] p;
        }
        /*note that, above procedure only execute in the local Leader process,
        *so we need below procedure to copy that to all other process
        */

        //broadcast the value of b to all the process in the same shooting
        MPI_Bcast(b, FfsBranch::size, MPI_INT, 0, FfsBranch::commLocal);
        total=0;
        for (int i = 0; i < FfsBranch::size; i += 1) {
            total+=b[i];
        }
        int allSize=lambdaGlobal.size();
        MPI_Bcast(&allSize, 1, MPI_INT, 0, FfsBranch::commLocal);
        //similar process as above
        int *p=new int[allSize];
        if (local->isLeader) {
            //convert the lambdaGlobal to the p
            for (int i=0;i<allSize;i++) {
                p[i]=lambdaGlobal[i];
            }
        }
        MPI_Bcast(p, allSize, MPI_INT, 0, FfsBranch::commLocal);
        lambdaGlobal.clear();
        for (int i=0;i<allSize;i++) {
            lambdaGlobal.push_back(p[i]);
        }
        delete[] p;
    }
    const std::string getName(int x) const {
        x%=total;
        int i;
        for (i = 0; i < FfsBranch::size; i += 1) {
            if (x<b[i]) {
                return generateName(x,i);
            }
            else {
                x-=b[i];
            }
        }
    }
    int getLambda(int x) const {
        return lambdaGlobal[x%total];
    }
    int getTotal() const {
        return total;
    }
private:
    int layer;

    //a and b is an array used for storing the size of lambda vector, a[i] is the size of lambda vector in the i shooting
    //a will only get part of the lambda value, which means, a only possess 1 shooting data, while b contains all shooting data, access by the branchId
    int *a,*b;
    //lambdaLocal only stores the lambda value in the certain layer and shooting
    std::vector<int> lambdaLocal;
    //lambdaGlobal stores the lambda value in all shooting in the certain layer
    std::vector<int> lambdaGlobal;
    int total;
    //generate a string with format "layer__branchId_inputparameter", and return it
    const std::string generateName(int x, int branchId=-1) const {
        if (branchId==-1) {
            branchId=local->id;
        }
        static char c[100];
        sprintf(c,"%d__%d_%d",layer,branchId,x);
        return c;
    }
};

//initialize the ffs process
bool ffsRequested(int argc, char **argv) {
    int i;
    for (i=0;i<argc;i++) {
        //search for the string "-ffs"
        if (strcmp(argv[i],"-ffs")==0) {
            //judge if there's 2 parameters after "-ffs", exemple: -ffs 512 ffs.input
            if (i+2<argc) {
                int numberShooting;
                //assign the first parameter to numberShooting
                if (sscanf(argv[i+1],"%d",&numberShooting)==1) {
                    MPI_Comm worldComm;
                    MPI_Comm_dup(MPI_COMM_WORLD,&worldComm);
                    //world is an instance which stores the communicator with each process is ffs
                    world=new MpiInfo(worldComm);
                    if (world->size>0&&numberShooting>0&&world->size%numberShooting==0) {
                        //the cpu number assigned to each shooting
                        int cpuEach=world->size/numberShooting;
                        //which shooting the current procedure need to handle
                        int myShooting=world->rank/cpuEach;
                        MPI_Comm localComm;
                        //seperate the process by myShooting
                        MPI_Comm_split(world->comm,myShooting,world->rank,&localComm);
                        local=new MpiInfo(localComm,myShooting);
                    }
                }
                else {
                    return false;
                }
                //read the file of ffs.input
                ffsParams=new FfsFileReader(argv[i+2]);
                return true;
            }
        }
    }
    return false;
}

//
class FfsRandomGenerator: public FfsBranch {
    public:
        //generate random seed in the leader process of local
        FfsRandomGenerator() {
            if (!local->isLeader) {
                return ;
            }
            std::srand(std::time(0) + world->rank * 1234567);
            std::rand();
            std::rand();
        }
        //generalize random number in main process and distrubute in all process
        int get() {
            int x;
            if (local->isLeader) {
                x=std::rand();
            }
            MPI_Bcast(&x, 1, MPI_INT, 0, FfsBranch::commLocal);
            return x;
        }
};

//set the velocity of atoms in gaussian distribution
int createVelocity(LAMMPS *lammps, const std::string &groupName, int temp, FfsRandomGenerator *pRng) {
    static char str[100];
    int seed=pRng->get();
    sprintf(str,"velocity %s create %d %d dist gaussian", groupName.c_str(), temp, seed);
    //set the velocity of the atoms
    lammps_command(lammps,str);
    return seed;
};
void runBatch(LAMMPS *lammps) {
    static char str[100];
    static bool inited=false;
    if (!inited) {
        //get the value of parameter "check_every"
        int everyStep=ffsParams->getInt("check_every");
        sprintf(str,"run %d pre no post no",everyStep);
        inited=true;
    }
    if (lammps) {
        //run check_every step
        lammps_command(lammps,str);
    }
}

void printBox(void *lmp, const std::string &xyzFinal) {
  if (local->isLeader) {
    static double low[3], high[3], xy, yz, xz;
    static int periodicity[3], boxChange;
    lammps_extract_box(lmp, low, high, &xy, &yz, &xz, periodicity, &boxChange);
    printf("%s (%f, %f, %f) - (%f, %f, %f)\n", xyzFinal.c_str(), low[0], low[1], low[2], high[0], high[1], high[2]);
  }
}

//print the status randomly according to the print_every
void printStatus(const int print_every, const int timestep, const int current, const int target) {
  if (local->isLeader) {
    if (std::rand() < 1.0 * RAND_MAX / print_every) {
      printf("[date=%d] [universe=%d] [steps=%d] : %d ... %d\n", std::time(0), local->id, timestep, current, target);
    }
  }
}

/**
*\param argc the number of parameters
*\param argv a pointer to the pointer of char, used for store the value of parameter
*/
int ffs_main(int argc, char **argv) {
    //initialization, 0 can be treated as nullptr
    runBatch(0);
    LAMMPS *lammps=new LAMMPS(argc,argv,local->comm);
    //process all the lammps input file
    lammps->input->file();
    //lammps_command(lammps,(char *)"set group all image 0 0 0");
    FfsTrajectoryReader continuedTrajectory;
    FfsTrajectoryWriter fileTrajectory;
    //get the value of the specific parameter
    int temperatureMean=ffsPar ams->getInt("temperature");
    const std::string waterGroupName = ffsParams->getString("water_group");
    int equilibriumSteps=ffsParams->getInt("equilibrium");
    int print_every = ffsParams->getInt("print_every");
    const std::vector<int> config_each_lambda = ffsParams->getVector("config_each_lambda");  
    const std::vector<int> lambdaList=ffsParams->getVector("lambda");
    static int lambda_A=lambdaList[0];
    FfsRandomGenerator rng;
    FfsFileTree *lastTree,*currentTree;

    //set velocity of the atoms and get the seed
	int velocitySeed=createVelocity(lammps, waterGroupName, temperatureMean, &rng);
    //parameter is number of configurations collected at each interface, if there's continue file, then continue the process
	FfsCountdown *fcd = new FfsCountdown(config_each_lambda[0] - continuedTrajectory.countPrecalculated(0)); 
	lammps_command(lammps,(char *)"run 0 pre yes post no");
	lastTree=0;
	currentTree=new FfsFileTree(&continuedTrajectory,0);
	while (1) {
		if (local->id == 0) {
		  sleep(1);
          //receive the file information from all other process
		  fileTrajectory.check();
		  if (!fcd->next()) {
			delete fcd;
			break;
		  }
		  continue;
		}
		bool ready=false;
		int lambda;
		while (1) {
            //run check_every steps
			runBatch(lammps);
			const double *lambdaReuslt=(const double *)lammps_extract_compute(lammps,(char *)"lambda",0,1);
			lambda=(int)lambdaReuslt[0];
			static int lambda_0=lambdaList[1];
            //update is a member variant with class "update" in lammps object, and ntimestep stores the timestep now 
			const int64_t timestep = lammps->update->ntimestep;
			printStatus(print_every, timestep, lambda, lambda_0);
			if (lambda<=lambda_A) {
				ready=true;
			}
            //when the lambda evolved to be smaller than lambdaA and larger than lambda0, break loop
			if (ready&&lambda>=lambda_0) {
				ready=false;
				break;
			}
			fileTrajectory.check();
            //if has reached the configuration number, then break loop
			if (!fcd->next()) {
				break;
			}
		}
		if (!fcd->next()) {
			delete fcd;
			break;
		}
		int64_t timestep=lammps->update->ntimestep;
        //if the timestep hasn't reach the equilibriumsteps, then rerun the loop
		if (timestep <= equilibriumSteps) {
		  continue;
		}
		static char strDump[100];
        //xyzFinal is a string that can reflect the layer number and current lambda value, 0__4_0 stands for NO.0 layer, branchId 4, and 0 stands for the size of lambda vector in the branch
		const std::string xyzFinal=currentTree->add(lambda);
		sprintf(strDump,"write_dump all xyz pool/xyz.%s",xyzFinal.c_str());
        /*
        *write the trajectory information into the file "trajectory.out.txt", exemple: "  ___ (__________)  >==1777855480               106360==>   40 (xyz.0__4_0)"
        *here the 1777855480 stands for velocity seed, and 106360 stands for the timestep number, 40 is the value of lambda
        */
		fileTrajectory.writeln((const char *)0,0,velocitySeed,timestep,xyzFinal.c_str(),lambda);
		lammps_command(lammps,strDump);
        //print the parameter of box
		printBox(lammps, xyzFinal);
		fcd->done();
	}
	lammps_command(lammps,(char *)"run 0 pre no post yes");
    
    //the second part, loop until finish
    const int n=lambdaList.size();
    for (int i=1;i+1<n;i++) {
        delete lastTree;
        lastTree=currentTree;
        currentTree=new FfsFileTree(&continuedTrajectory,i);
        lastTree->commit();
        FfsCountdown *fcd = new FfsCountdown(config_each_lambda[i] - continuedTrajectory.countPrecalculated(i));  
        const int lambda_next=lambdaList[i+1];
        while (1) {
            //the process with id 0 only possess file information
            if (local->id == 0) {
              sleep(0.5);
              fileTrajectory.check();
              if (!fcd->next()) {
                delete fcd;
                break;
              }
              continue;
            }
            static char strReadData[100];
            const int initConfig=rng.get();
            const std::string xyzInit=lastTree->getName(initConfig);
            const int lambdaInit=lastTree->getLambda(initConfig);
            //get the configuration from stored data
            sprintf(strReadData,"read_dump pool/xyz.%s 0 x y z box no format xyz",xyzInit.c_str());
            lammps_command(lammps,strReadData);
            //regenerate seeds
            int velocitySeed=createVelocity(lammps, waterGroupName, temperatureMean, &rng);
            if (local->isLeader) {
                printf("[date=%d] [universe=%d] [initialFile=%s] [velocitySeed=%d]\n", std::time(0), local->id, xyzInit.c_str(), velocitySeed);
            }
            lammps_command(lammps,(char *)"run 0 pre yes post no");
            int lambda_calc;
            while (1) {
                runBatch(lammps);
                const double *lambdaReuslt=(const double *)lammps_extract_compute(lammps,(char *)"lambda",0,1);
                lambda_calc=(int)lambdaReuslt[0];
                const int64_t timestep = lammps->update->ntimestep;
                printStatus(print_every, timestep, lambda_calc, lambda_next);
                if (lambda_calc<=lambda_A||lambda_calc>=lambda_next) {
                    break;
                }
                fileTrajectory.check();
                if (!fcd->next()) {
                    break;
                }
            }
            lammps_command(lammps,(char *)"run 0 pre no post yes");
            if (!fcd->next()) {
                delete fcd;
                break;
            }
            int64_t timestep=lammps->update->ntimestep;
            //the system has returned to the initial point, so print and continue
            if (lambda_calc<=lambda_A) {
              if (local->isLeader) {
                printf("%3d (xyz.%s)  >==%010d %20lld==>  %3d (__________)\n", lambdaInit, xyzInit.c_str(), velocitySeed, timestep, lambda_calc);
              }
              continue;
            }
            //if reach the next layer, store it
            if (lambda_calc>=lambda_next) {
                static char strDump[100];
                const std::string xyzFinal=currentTree->add(lambda_calc);
                sprintf(strDump,"write_dump all xyz pool/xyz.%s",xyzFinal.c_str());
                fileTrajectory.writeln(xyzInit.c_str(),lambdaInit,velocitySeed,timestep,xyzFinal.c_str(),lambda_calc);
                lammps_command(lammps,strDump);
                printBox(lammps, xyzFinal);
                fcd->done();
                continue;
            }
        }
    }
    delete lammps;
    delete local;
    delete world;
    return 0;
}

