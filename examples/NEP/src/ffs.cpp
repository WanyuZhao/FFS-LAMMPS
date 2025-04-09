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
struct MpiInfo {
    const int LEADER;
    const int id;
    MPI_Comm comm;
    int size,rank;
    bool isLeader;
    MpiInfo(MPI_Comm comm,int id=0):comm(comm),LEADER(0),id(id) {
        MPI_Comm_size(comm,&size);
        MPI_Comm_rank(comm,&rank);
        isLeader=LEADER==rank;
    };
};
const MpiInfo *world,*local;
class FfsBranch {
public:
    FfsBranch() {
      if (!FfsBranch::commInited) {
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
MPI_Comm FfsBranch::commLeader,FfsBranch::commLocal;
int FfsBranch::size=0;
struct FfsFileReader: public FfsBranch {
    std::map<std::string,std::string> dict;
    FfsFileReader(const char *filename) {
        if (world->isLeader) {
            FILE *f=fopen(filename,"r");
            static char cLine[1024];
            while (fgets(cLine,sizeof(cLine),f)) {
                std::string sLine(cLine);
                int l=sLine.find_first_of('#');
                if (l!=std::string::npos) {
                    sLine=sLine.substr(0,l);
                }
                const std::string SPACE=" \t\n\v\f\r";
                int begin1=sLine.find_first_not_of(SPACE);
                if (begin1==std::string::npos) {
                    continue;
                }
                int end1=sLine.find_first_of(SPACE,begin1);
                if (end1==std::string::npos) {
                    end1=sLine.length();
                }
                std::string key=sLine.substr(begin1,end1-begin1);
                dict[key]="";
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
    int getInt(const std::string &name) const {
        int y=0;
        const std::string v=getString(name);
        sscanf(v.c_str(),"%d",&y);
        return y;
    }
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
    const std::string getString(const std::string &name) const {
        char *p;
        int l;
        if (world->isLeader) {
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
    void writeln(const char *format,...) {
        if (!local->isLeader) {
            return ;
        }
        static char buffer[MAX_LENGTH+1];
        va_list args;
        va_start(args,format);
        vsprintf(buffer,format,args);
        va_end(args);
        putstr(buffer);
    }
    void check() {
        if (!world->isLeader) {
            return ;
        }
        static char buffer[MAX_LENGTH];
        while (1) {
            int flag;
            MPI_Status status;
            MPI_Iprobe(MPI_ANY_SOURCE, FfsBranch::TAG_FILEWRITER_LINE, FfsBranch::commLeader, &flag, &status);
            if (flag) {
                int l;
                MPI_Get_count(&status,MPI_CHAR,&l);
                printf("[date=%d] world leader will receive TAG_FILEWRITER_LINE from %d\n", std::time(0), status.MPI_SOURCE);
                MPI_Recv(buffer, l, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, FfsBranch::commLeader, &status);
                printf("[date=%d] world leader did receive TAG_FILEWRITER_LINE from %d\n", std::time(0), status.MPI_SOURCE);
                putstr0(status.MPI_SOURCE,buffer);
            }
            else {
                break;
            }
        }
    }
    void pushGroup(int n) {
        nFlush+=n;
    }
private:
    FILE *f;
    void putstr(char *s) {
        s[MAX_LENGTH-1]='\0';
        if (world->isLeader) {
            putstr0(0,s);
        }
        else {
            int l=strlen(s);
            printf("[date=%d] universe %d will send TAG_FILEWRITER_LINE\n", std::time(0), local->id);
            MPI_Send(s, l + 1, MPI_CHAR, 0, FfsBranch::TAG_FILEWRITER_LINE, FfsBranch::commLeader);
            printf("[date=%d] universe %d did send TAG_FILEWRITER_LINE\n", std::time(0), local->id);
        }
    }
    int nFlush;
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
    void writeln(const char *xyzInit, const char *xyzFinal,const std::vector<int> &v) {
        int n=v.size();
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
            while (1) {
                int lambda,layer,branch,count;
                int ret=fscanf(f,"%*s%*s%*s%*s%*s%d%*[^.]%*c%d%*c%*c%d%*c%d%*[^\n]%*c",&lambda,&layer,&branch,&count);
                if (ret==EOF) {
                    break;
                }
                v[branch].push_back(layer);
                v[branch].push_back(count);
                v[branch].push_back(lambda);
            }
            fclose(f);
            int i,j;
            if (1) {
                const std::vector<int> vv=v[0];
                n=vv.size();
                p=new int[n];
                std::copy(vv.begin(),vv.end(),p);
            }
            for (int i=1;i<(int)v.size();i++) {
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
            MPI_Probe(0, FfsBranch::TAG_FILEREADER, FfsBranch::commLeader, &status);
            MPI_Get_count(&status,MPI_INT,&n);
            p=new int[n];
            printf("[date=%d] universe %d will receive TAG_FILEREADER\n", std::time(0), local->id);
            MPI_Recv(p, n, MPI_INT, 0, FfsBranch::TAG_FILEREADER, FfsBranch::commLeader, &status);
            printf("[date=%d] universe %d did receive TAG_FILEREADER\n", std::time(0), local->id);
        }
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
        MPI_Barrier(FfsBranch::commLocal);
    }
    const std::vector<int> &get(int layer) const {
        if (layer<lambdaLocal.size()) {
            return lambdaLocal[layer];
        }
        else {
            return emptyVector;
        }
    }
    int countPrecalculated(int layer) const {
        int x;
        if (local->isLeader) {
            int t=get(layer).size();
            MPI_Allreduce(&t, &x, 1, MPI_INT, MPI_SUM, FfsBranch::commLeader);
        }
        MPI_Bcast(&x, 1, MPI_INT, 0, FfsBranch::commLocal);
        return x;
    }
private:
    std::vector< std::vector<int> > lambdaLocal;
    std::vector<int> emptyVector;
};
class FfsCountdown: public FfsBranch {
public:
    FfsCountdown(int n) {
        remains=n;
        terminated=n<=0;
    }
    ~FfsCountdown() {
    }
    void done(int n=1) {
        if (!local->isLeader) {
            return ;
        }
        if (terminated) {
            return ;
        }
        int x=n;
        if (world->isLeader) {
            remains-=x;
        }
        else {
            printf("[date=%d] universe %d will send TAG_COUNTDOWN_DONE\n", std::time(0), local->id);
            MPI_Send(&x, 1, MPI_INT, 0, FfsBranch::TAG_COUNTDOWN_DONE, FfsBranch::commLeader);
            printf("[date=%d] universe %d did send TAG_COUNTDOWN_DONE\n", std::time(0), local->id);
        }
    }
    bool next() {
        if (terminated) {
            return false;
        }
        int ret;
        if (local->isLeader) {
            if (world->isLeader) {
                while (1) {
                    int flag;
                    MPI_Status status;
                    MPI_Iprobe(MPI_ANY_SOURCE, FfsBranch::TAG_COUNTDOWN_DONE, FfsBranch::commLeader, &flag, &status);
                    if (flag) {
                        int x;
                        printf("[date=%d] world leader will receive TAG_COUNTDOWN_DONE from %d\n", std::time(0), status.MPI_SOURCE);
                        MPI_Recv(&x, 1, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, FfsBranch::commLeader, &status);
                        printf("[date=%d] world leader did receive TAG_COUNTDOWN_DONE from %d\n", std::time(0), status.MPI_SOURCE);
                        remains-=x;
                    }
                    else {
                        break;
                    }
                }
                if (remains<=0) {
                    terminated=true;
                    int i;
                    for (i = 1; i < FfsBranch::size; i += 1) {
                        printf("[date=%d] world leader will send TAG_COUNTDOWN_TERMINATE to %d\n", std::time(0), i);
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
        MPI_Bcast(&ret, 1, MPI_INT, 0, FfsBranch::commLocal);
        if (ret==0) {
            terminated=true;
        }
        return ret;
    }
private:
    bool terminated;
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
            lambdaLocal=ftr->get(layer);
            a[local->id]=lambdaLocal.size();
        }
    }
    ~FfsFileTree() {
        delete[] b;
        delete[] a;
    }
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
    void commit() {
        if (local->isLeader) {
            MPI_Allreduce(a, b, FfsBranch::size, MPI_INT, MPI_SUM, FfsBranch::commLeader);
            int mySize=lambdaLocal.size();
            int currentSize;
            int maxSize;
            MPI_Allreduce(&mySize, &maxSize, 1, MPI_INT, MPI_MAX, FfsBranch::commLeader);
            int *p=new int[maxSize];
            int i,j;
            for (i = 0; i < FfsBranch::size; i += 1) {
                if (i==local->id) {
                    currentSize=mySize;
                    std::copy(lambdaLocal.begin(),lambdaLocal.end(),p);
                }
                MPI_Bcast(&currentSize, 1, MPI_INT, i, FfsBranch::commLeader);
                MPI_Bcast(p, currentSize, MPI_INT, i, FfsBranch::commLeader);
                for (j=0;j<currentSize;j++) {
                    lambdaGlobal.push_back(p[j]);
                }
            }
            delete[] p;
        }
        MPI_Bcast(b, FfsBranch::size, MPI_INT, 0, FfsBranch::commLocal);
        total=0;
        for (int i = 0; i < FfsBranch::size; i += 1) {
            total+=b[i];
        }
        int allSize=lambdaGlobal.size();
        MPI_Bcast(&allSize, 1, MPI_INT, 0, FfsBranch::commLocal);
        int *p=new int[allSize];
        if (local->isLeader) {
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
    int *a,*b;
    std::vector<int> lambdaLocal;
    std::vector<int> lambdaGlobal;
    int total;
    const std::string generateName(int x, int branchId=-1) const {
        if (branchId==-1) {
            branchId=local->id;
        }
        static char c[100];
        sprintf(c,"%d__%d_%d",layer,branchId,x);
        return c;
    }
};
bool ffsRequested(int argc, char **argv) {
    int i;
    for (i=0;i<argc;i++) {
        if (strcmp(argv[i],"-ffs")==0) {
            if (i+2<argc) {
                int numberShooting;
                if (sscanf(argv[i+1],"%d",&numberShooting)==1) {
                    MPI_Comm worldComm;
                    MPI_Comm_dup(MPI_COMM_WORLD,&worldComm);
                    world=new MpiInfo(worldComm);
                    if (world->size>0&&numberShooting>0&&world->size%numberShooting==0) {
                        int cpuEach=world->size/numberShooting;
                        int myShooting=world->rank/cpuEach;
                        MPI_Comm localComm;
                        MPI_Comm_split(world->comm,myShooting,world->rank,&localComm);
                        local=new MpiInfo(localComm,myShooting);
                    }
                }
                else {
                    return false;
                }
                ffsParams=new FfsFileReader(argv[i+2]);
                return true;
            }
        }
    }
    return false;
}
class FfsRandomGenerator: public FfsBranch {
    public:
        FfsRandomGenerator() {
            if (!local->isLeader) {
                return ;
            }
            std::srand(std::time(0) + world->rank * 1234567);
            std::rand();
            std::rand();
        }
        int get() {
            int x;
            if (local->isLeader) {
                x=std::rand();
            }
            MPI_Bcast(&x, 1, MPI_INT, 0, FfsBranch::commLocal);
            return x;
        }
};

int createVelocity(LAMMPS *lammps, const std::string &groupName, int temp, FfsRandomGenerator *pRng) {
    static char str[100];
    int seed=pRng->get();
    sprintf(str,"velocity %s create %d %d dist gaussian", groupName.c_str(), temp, seed);
    lammps_command(lammps,str);
    return seed;
};
void runBatch(LAMMPS *lammps) {
    static char str[100];
    static bool inited=false;
    if (!inited) {
        int everyStep=ffsParams->getInt("check_every");
        sprintf(str,"run %d pre no post no",everyStep);
        inited=true;
    }
    if (lammps) {
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

void printStatus(const int print_every, const int timestep, const int current, const int target) {
  if (local->isLeader) {
    if (std::rand() < 1.0 * RAND_MAX / print_every) {
      printf("[date=%d] [universe=%d] [steps=%d] : %d ... %d\n", std::time(0), local->id, timestep, current, target);
    }
  }
}

int ffs_main(int argc, char **argv) {
    runBatch(0);
    LAMMPS *lammps=new LAMMPS(argc,argv,local->comm);
    lammps->input->file();
    //lammps_command(lammps,(char *)"set group all image 0 0 0");
    FfsTrajectoryReader continuedTrajectory;
    FfsTrajectoryWriter fileTrajectory;
    int temperatureMean=ffsParams->getInt("temperature");
    const std::string waterGroupName = ffsParams->getString("water_group");
    int equilibriumSteps=ffsParams->getInt("equilibrium");
    int print_every = ffsParams->getInt("print_every");
    const std::vector<int> config_each_lambda = ffsParams->getVector("config_each_lambda");   
    const std::vector<int> lambdaList=ffsParams->getVector("lambda");
    static int lambda_A=lambdaList[0];
    FfsRandomGenerator rng;
    FfsFileTree *lastTree,*currentTree;

    int velocitySeed=createVelocity(lammps, waterGroupName, temperatureMean, &rng);
    FfsCountdown *fcd = new FfsCountdown(config_each_lambda[0] - continuedTrajectory.countPrecalculated(0));  
    lammps_command(lammps,(char *)"run 0 pre yes post no");
    lastTree=0;
    currentTree=new FfsFileTree(&continuedTrajectory,0);
    while (1) {
        if (local->id == 0) {
          sleep(1);
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
            runBatch(lammps);
            const double *lambdaReuslt=(const double *)lammps_extract_compute(lammps,(char *)"lambda",0,1);
            lambda=(int)lambdaReuslt[0];
            static int lambda_0=lambdaList[1];
            const int64_t timestep = lammps->update->ntimestep;
            printStatus(print_every, timestep, lambda, lambda_0);
            if (lambda<=lambda_A) {
                ready=true;
            }
            if (ready&&lambda>=lambda_0) {
                ready=false;
                break;
            }
            fileTrajectory.check();
            if (!fcd->next()) {
                break;
            }
        }
        if (!fcd->next()) {
            delete fcd;
            break;
        }
        int64_t timestep=lammps->update->ntimestep;
        if (timestep <= equilibriumSteps) {
          continue;
        }
        static char strDump[100];
        const std::string xyzFinal=currentTree->add(lambda);
        sprintf(strDump,"write_dump all xyz pool/xyz.%s",xyzFinal.c_str());
        fileTrajectory.writeln((const char *)0,0,velocitySeed,timestep,xyzFinal.c_str(),lambda);
        lammps_command(lammps,strDump);
        printBox(lammps, xyzFinal);
        fcd->done();
    }
    lammps_command(lammps,(char *)"run 0 pre no post yes");


    const int n=lambdaList.size();
    for (int i=1;i+1<n;i++) {
        delete lastTree;
        lastTree=currentTree;
        currentTree=new FfsFileTree(&continuedTrajectory,i);
        lastTree->commit();
        FfsCountdown *fcd = new FfsCountdown(config_each_lambda[i] - continuedTrajectory.countPrecalculated(i));   
        const int lambda_next=lambdaList[i+1];
        while (1) {
            if (local->id == 0) {
              sleep(1);
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
            sprintf(strReadData,"read_dump pool/xyz.%s 0 x y z box no format xyz",xyzInit.c_str());
            lammps_command(lammps,strReadData);
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
            if (lambda_calc<=lambda_A) {
              if (local->isLeader) {
                printf("%3d (xyz.%s)  >==%010d %20lld==>  %3d (__________)\n", lambdaInit, xyzInit.c_str(), velocitySeed, timestep, lambda_calc);
              }
              continue;
            }
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

