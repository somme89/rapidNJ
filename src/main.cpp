#include "stdinclude.h"
#include "rdDataInitialiser.h"
#include "distMatrixReader.hpp"
#include "cluster_pair.h"
#include "rapidNJ.h"
#include "rapidNJDisk.h"
#include "rapidNJMem.hpp"
#include "simpleNJ.h"
#include "dataloader.hpp"
#include "dataloaderStockholm.hpp"
#include "dataloaderFasta.hpp"
#include "dataloaderPhylip.hpp"
#include "JCdistance.hpp"
#include "hammingDistance.hpp"
#include "KimuraDistance.hpp"
#include "getopt_pp/getopt_pp.h"
#include "bitStringUtils.hpp"
#include <iomanip>
#include <math.h>

using namespace std;

const string VERSION = "2.2.3";

#ifdef __APPLE__
#include <sys/types.h>
#include <sys/sysctl.h>
#endif

#ifdef __linux__
#include <sys/time.h>
timeval startTime,endTime;
#endif

#ifdef __WINDOWS__
#include <Windows.h>
#endif

bool distanceMatrixInput = true;
int matrixSize = -1;
int numCores = 1;

namespace options {
  bool verbose;
  string fileName;
  int memSize;
  int cores;
  string cacheDir;
  string percentageMemoryUsage;
  string distMethod;
  string inputFormat;
  string outputFormat;
  bool fastdist;
  int replicates;
  string inputtype;
  bool rapidNJ;
  bool simpleNJ;
  bool gpu;
  bool negative_branches;
  string outputFile;  
};


class distMatrixData {
public:
  distType** matrix;
  vector<string>* sequenceNames;
  diskMatrix* dm;
};

#ifdef __linux__

void printTime(string message){
  double startm = (startTime.tv_sec*1000.0 + startTime.tv_usec/1000.0)/1000;
  gettimeofday(&endTime,NULL);
  double endm = (endTime.tv_sec*1000.0 + endTime.tv_usec/1000.0)/1000;
  cerr << message << " " << (endm -startm ) << " sec \n";
}

void printTime(){
  double startm = (startTime.tv_sec*1000.0 + startTime.tv_usec/1000.0)/1000;
  gettimeofday(&endTime,NULL);
  double endm = (endTime.tv_sec*1000.0 + endTime.tv_usec/1000.0)/1000;
  cerr << (endm -startm );
}
#endif

int getMatrixSize(){
  ifstream is;
  is.open(options::fileName.data(),ifstream::in);  
  if(!is.is_open()){
    cerr << "Could not read file: " << options::fileName << "\n";
    exit(1);
  }
  char buf[256];
  is.read(buf, 256);
  is.close();
  return atoi(buf);
}

string guessPhylipType(string filename){
  ifstream is;
  is.open(filename.data(),ifstream::in);  
  if(!is.is_open()){
    cerr << "Could not read file: " << filename << "\n";
    exit(1);
  }
  char buf[256];
  is.read(buf, 256);
  is.close();
  int i = 0;
  bool foundNumber = false;
  for(; i < 256; i++){
    if(!foundNumber && buf[i] > 32) {
      foundNumber = true;
    } else if(foundNumber && (buf[i] == 32 || buf[i] == 11)){
      break;
    } else if(buf[i] == 10){
      //phylip distance matrix
      return "pd";
    }
  }  
  for(; i < 256; i++){
    if(buf[i] > 32){
      //alignment
      return "pa";
    } else if(buf[i] == 10){
      //distance matrix
      return "pd";
    }
  }
  return "pd";
}


void configureNumberOfCores(){
  // Configure number of cores to use
#ifdef __WINDOWS__
  SYSTEM_INFO sysinfo;
  GetSystemInfo( &sysinfo );
  numCores = sysinfo.dwNumberOfProcessors;  
#endif
#if defined(HW_NCPU) || defined(HW_AVAILCPU) 
  int mib[4];
  size_t len = sizeof(numCores);

  /* set the mib for hw.ncpu */
  mib[0] = CTL_HW;
  mib[1] = HW_AVAILCPU;  

  /* get the number of CPUs from the system */
  sysctl(mib, 2, &numCores, &len, NULL, 0);

  if( numCores < 1 ) {
    //try HW_NCPU;
    mib[1] = HW_NCPU;
    sysctl( mib, 2, &numCores, &len, NULL, 0 );
  }
#endif
#ifdef _SC_NPROCESSORS_ONLN
  numCores= sysconf(_SC_NPROCESSORS_ONLN);
#endif
  if(options::cores > 0){
    if(options::cores > numCores){
      cerr << "WARNING: the system only have " << numCores << ". Using " << options::cores << " cores will result in suboptimal performance." << endl;
    }
    numCores = options::cores;
  }
  if(numCores < 1) {
    numCores = 1;
  }
  if (options::verbose) {
    cerr << "Using " << numCores << " core(s) for distance estimation" << endl;
  }
}

/*Returns the system memory size in MB */
double getMemSize(){
  double retval = -1;
  if(options::memSize != -1){
    retval = options::memSize * 1024.0 * 1024.0;
  } else {
#ifdef __APPLE__
  int mib[2];
  int64_t physical_memory;
  size_t length;
  mib[0] = CTL_HW;
  mib[1] = HW_MEMSIZE;
  length = sizeof(int64_t);
  sysctl(mib, 2, &physical_memory, &length, NULL, 0);
  retval = (double) physical_memory;
#elif defined __WINDOWS__
  MEMORYSTATUSEX status;
  status.dwLength = sizeof(status);
  GlobalMemoryStatusEx(&status);
  retval = (double) status.ullTotalPhys;
#else
  double pages = (double) sysconf(_SC_PHYS_PAGES);
  double page_size = (double) sysconf(_SC_PAGE_SIZE);
  retval = pages * page_size;  
#endif
  }
  if(retval > 0.0) {
    return retval * 0.8;
  } else {
    cerr << "Could not read memory size. Please state the amount of available memory using the '-m'option." << endl;
    exit(1);
    return -1;
  }
} 

void printDistanceMatrix(ostream& out, distMatrixData* data) {
  out << "\t" << matrixSize << endl;
  for (int i = 0; i < matrixSize; i++) {
    out << (*data->sequenceNames)[i] << "\t";
    for (int j = 0; j < matrixSize; j++) {
      out << setprecision(6) << fixed << data->matrix[i][j] << " ";
    }
    out << endl;
  }
}

void printDistanceMatrixDisk(ostream& out, distMatrixData* data) {
  out << "\t" << matrixSize << endl;
  distType* row = new distType[matrixSize];
  for (int i = 0; i < matrixSize; i++) {
    out << (*data->sequenceNames)[i] << "\t";
    for (int j = 0; j < matrixSize; j++) {
      data->dm->readArray(row, i, matrixSize);
      out << setprecision(6) << fixed << row[j] << " ";
    }
    out << endl;
  }
  delete[] row;
}

distMatrixData* computeDistanceMatrix(bool useDiskMatrix, ostream &out, bool printMatrix, dataloader* dl) {  
  /*if(options::gpu){
  computeDistanceMatrixGPU();
  return;
  }*/  

  if(options::fastdist && options::verbose){
    cerr << "Fastdist is enabled" << endl;
  }

  distMatrixData* retVal = new distMatrixData();  
  if(useDiskMatrix) {
    retVal->dm = new diskMatrix(options::cacheDir, matrixSize);
  }

  retVal->sequenceNames = dl->getSequenceNames();
  matrixSize = dl->getSequenceCount();
  // process data
  if(options::distMethod == "jc"){    
    if (options::verbose) {
      cerr << "Using JC algorithm to calculate distances" << endl;
    }
    JCdistance* alg = new JCdistance(options::verbose, options::fastdist, dl, retVal->dm);
    alg->computeDistanceMatrix(numCores);
    retVal->matrix = alg->getDistanceMatrix();
    delete alg;
  } else if(options::distMethod == "kim" || options::distMethod == ""){
    if (options::verbose) {
      cerr << "Using Kimura algorithm to calculate distances" << endl;
    }
    KimuraDistance* alg = new KimuraDistance(options::verbose, options::fastdist, dl, retVal->dm);    
    alg->computeDistances(numCores);
    retVal->matrix = alg->getDistanceMatrix();
    delete alg;
  } else {
    cerr << "ERROR: Unknown sequence evolution model" << endl;
    exit(1);
  }

  if(printMatrix) {
    //Print the matrix
    if(retVal->dm == NULL) {
      printDistanceMatrix(out, retVal);
    } else {      
      printDistanceMatrixDisk(out, retVal);
      delete retVal->dm;
    }
  }
  return retVal;
}

polytree* runRapidNJ(int matrixSize, distMatrixReader* reader, ProgressBar* pb){   
  if(options::verbose){
    cerr << "Computing phylogetic tree... \n";
  }
  rapidNJ* sorted = new rapidNJ(reader, matrixSize, options::negative_branches, pb);
  polytree* tree = sorted->run();
  //TODO delete rd
  delete sorted;
  return tree;
}

polytree* runSimpleNJ(distMatrixReader* reader, ProgressBar* pb){
  simpleNJ* njs = new simpleNJ(reader, matrixSize, options::negative_branches, pb);
  polytree* tree = njs->run();
  delete njs;
  return tree;
}

polytree* runRapidMNJ(int sortedMatrixSize, distMatrixReader* reader, ProgressBar* pb){
  if(options::verbose){
    cerr << "Computing phylogetic tree... \n";
  }  
  rapidNJMem* nj = new rapidNJMem(reader,matrixSize,sortedMatrixSize,options::verbose, options::negative_branches, pb);
  polytree* tree = nj->run();
  delete nj;
  return tree;
}

polytree* runDiskNJ(ostream &out, int datastructureSize, dataloader* dl, ProgressBar* pb){
  if(options::verbose){
    cerr << "Reading data... \n";
  }
  rdDataInitialiser* reader;

  if(distanceMatrixInput){
    reader = new rdDataInitialiser(options::verbose, datastructureSize, options::cacheDir, options::fileName);
    bool status = reader->read_data();
    if(!status){
      cerr << "Could not read distance matrix in file " << options::fileName << endl;
      exit(1);
    }
  } else {    
    distMatrixData* matrixData = computeDistanceMatrix(true, out, false, dl);    
    reader = new rdDataInitialiser(options::verbose, datastructureSize, options::cacheDir, matrixSize);
    reader->initializeFromExistingMatrix(matrixData->sequenceNames, matrixData->dm);    
    delete matrixData;
  }
   
  if(options::verbose){
    cerr << "Computing phylogetic tree... \n";
  }
  rapidNJDisk *rd = new rapidNJDisk(reader, options::verbose, options::negative_branches, pb);
  polytree* tree = rd->run();
  delete rd;  
  delete reader;
  return tree;
}

/*void computeDistanceMatrixGPU(){
dataloader* dl;
#ifndef ENABLEGPU
cerr << "GPU distance estimation is not enabled. Please recompile rapidNJ with CUDA support." << endl;
exit(1);
#endif
if(options::verbose){
cerr << "Using GPU to build distance matrix" << endl;
}
// load data 
if (options::verbose) {
cerr << "Reading data...\n";
}
InputType type = UNKNOWN;
if(options::inputtype == "d") type = DNA;
if(options::inputtype == "p") type = PROTEIN;
if(options::inputFormat == "sth"){
dl = new dataloaderStockholm();    
dl->readData(options::fileName, options::fastdist | options::gpu, options::verbose, true, type, true);
} else if(options::inputFormat == "fa") {
cerr << "GPU distance estimation with FASTA formatted input are not supported yet." << endl;
exit(1);
} else if(options::inputFormat == "pa") {
cerr << "GPU distance estimation with phylip formatted input are not supported yet." << endl;
exit(1);
} else {
cerr << "Unkown input format " << options::inputFormat << endl;
exit(0);
}
sequenceNames = dl->getSequenceNames();
matrixSize = dl->getSequenceCount();

// process data
if(options::distMethod == "kim" || options::distMethod == ""){
if (options::verbose) {
cerr << "Using Kimura algorithm to calculate distances" << endl;
}
KimuraDistance* alg = new KimuraDistance(options::verbose, true, dl,false);
alg->computeDistancesGPU();
distanceMatrix = alg->getDistanceMatrix();
} else {
cerr << "unknown method" << endl;
}
}*/

dataloader* loadAlignment(){
  // load data
  dataloader* dl;
  if (options::verbose) {
    cerr << "Reading data...\n";
  }
  InputType type = UNKNOWN;
  if(options::inputtype == "d") type = DNA;
  if(options::inputtype == "p") type = PROTEIN;
  if(options::inputFormat == "sth"){
    dl = new dataloaderStockholm();
    dl->readData(options::fileName, options::fastdist, options::verbose, type,false);
  } else if (options::inputFormat == "fa") {
    dl = new dataloaderFasta();
    dl->readData(options::fileName, options::fastdist, options::verbose, type,false);
  } else if(options::inputFormat == "pa"){
    //TODO fix this
    cerr << "Phylip alignments cannot be read yet. Use stockholm format for alignments."<< endl;
    exit(1);
    dl = new dataloaderPhylip();
    dl->readData(options::fileName, options::fastdist, options::verbose, type,false);
  } else {
    cerr << "Unkown input format " << options::inputFormat << endl;
    exit(0);
  } 
  return dl;
}

distMatrixReader* getDistanceMatrixData(ostream &out, bool halfMatrix, dataloader* dl) {  
  distMatrixReader* reader;
  if(distanceMatrixInput){
    reader = new distMatrixReader(options::verbose, options::fileName, matrixSize, halfMatrix);
    if(options::verbose){
      cerr << "Reading distance matrix... \n";
    }
    reader->read_data(NULL);
  } else {
    if(options::verbose){
      cerr << "Computing distance matrix... \n";
    }
    distMatrixData* matrixData = computeDistanceMatrix(false, out, false, dl);
    reader = new distMatrixReader(options::verbose, matrixSize, halfMatrix, matrixData->sequenceNames, matrixData->matrix);
    reader->initializeData();
    delete matrixData;
  }
  return reader;
}

polytree* computeTree(ostream &out, dataloader* dl, ProgressBar* pb){
  // Compute the memory requirements for the three different algorithms.
  bool autoDecide = true;
  double matrixSized = (double) matrixSize;
  double systemMemory = getMemSize();    
  double matrixMemUsage = ((double) sizeof(distType)) *  matrixSized * matrixSized; 
  double sortedMatrixMemUsage = matrixSized * matrixSized * ((double)sizeof(cluster_pair));
  int sortedMatrixSize = (int) ((systemMemory - (matrixMemUsage/2.0)) / (matrixSized*sizeof(cluster_pair)));
  sortedMatrixSize = min(sortedMatrixSize,matrixSize);
  int diskSortedMatrixSize = (int) (systemMemory / (distType)(matrixSize*(sizeof(cluster_pair)+sizeof(distType))));
  diskSortedMatrixSize = min(diskSortedMatrixSize,matrixSize);
  diskSortedMatrixSize = max(diskSortedMatrixSize, min(5, matrixSize));
  polytree* tree;

  // select which algorithm to use based on either parameters or memory requirements
  if(options::rapidNJ || options::cacheDir != "" || options::percentageMemoryUsage != "" || options::simpleNJ) {
    autoDecide = false;
  }
  if(options::verbose){
    cerr << "Matrix size: " << matrixSize << endl;
    cerr << (systemMemory / 1024 / 1024 / 0.8) << " MB of memory is available" << endl;
  }

  if(options::rapidNJ || (autoDecide && sortedMatrixMemUsage + matrixMemUsage <= systemMemory)) {
    if(options::verbose){
      cerr << "Using RapidNJ \n";
      cerr << "Using " << (matrixMemUsage /1024 /1024) << " MB for distance matrix" << endl;
      cerr << "Using " << (sortedMatrixMemUsage/1024/1024) << " MB for sortedMatrix" << endl;
      cerr << "Total memory consumption is " << (matrixMemUsage + sortedMatrixMemUsage)/1024/1024 << " MB" << endl;
    }
    if(sortedMatrixMemUsage > systemMemory - matrixMemUsage){
      cerr << "WARNING: There's not enough memory to use RapidNJ. Consider using another algorithm." << endl;
    }
    distMatrixReader* reader = getDistanceMatrixData(out, false, dl);
    tree = runRapidNJ(matrixSize, reader, pb);
  } else if(options::percentageMemoryUsage != ""  || (autoDecide && sortedMatrixSize >= matrixSize * MIN_SORTED_MATRIX_SIZE)){
    if(options::verbose){
      cerr << "Using Memory efficient RapidNJ \n";
    }
    if(options::percentageMemoryUsage != ""){
      // try to use the user supplied argument
      int percentage;
      percentage = atoi(options::percentageMemoryUsage.data());
      if(percentage < 0 || percentage > 100){
        cerr << "The memory use percentage must be >=0 and <=100 " << endl;
        exit(1);
      }
      int tempSize = (int) (matrixSize * (percentage / 100.0));
      if(tempSize > sortedMatrixSize){
        cerr << "WARNING: Not enough memory for " << percentage << "% of the sorted matrix. Reduce the size of the sorted matrix or use RapidDiskNJ." << endl;
      }
      sortedMatrixSize = tempSize;
    }
    if(sortedMatrixSize < matrixSize * MIN_SORTED_MATRIX_SIZE) {
      cerr << "WARNING: the amount of available memory is too low for the memory efficient RapidNJ algorithm to run efficiently. Consider using RapidDiskNJ." << endl;
    }
    if (sortedMatrixSize < 1) {
      sortedMatrixSize = 1;
    }
    if(options::verbose){     
      cerr << "Sorted matrix has " << sortedMatrixSize << " columns" << endl;
    }
    distMatrixReader* reader = getDistanceMatrixData(out, true, dl);
    tree = runRapidMNJ(sortedMatrixSize, reader, pb);
  } else if(options::simpleNJ) {
    if(options::verbose){
      cerr << "Using naive NJ \n";
    }   
    distMatrixReader* reader = getDistanceMatrixData(out, false, dl);
    tree = runSimpleNJ(reader, pb);
  } else {
    if(options::verbose) {
      cerr << "Using RapidDiskNJ algorithm\n";
      cerr << "Sorted matrix has " << diskSortedMatrixSize << " columns" << endl;
    }    
    tree = runDiskNJ(out, diskSortedMatrixSize, dl, pb);
  }
  return tree;  
}

void bootstrapTree(ostream &out, polytree* tree, dataloader* dl, ProgressBar* pb) {
  for(int i = 0; i < options::replicates; i++) {
    dl->sample_sequences();
    pb->childProgress(1.0 / (options::replicates + 1.0));
    polytree* replicate = computeTree(out, dl, pb);
    if(options::verbose) {
      cerr << "Comparing trees..." << endl;
    }
	//cout << "---------------------" << i << "-------------------------" << endl;
    tree->compareTreeBootstrap(replicate);
    delete replicate;
  }
  cout << endl;
}

void printUsage(){
  cerr << "Rapid neighbour-joining. An implementation of the canonical neighbour-joining method which utilize a fast search heuristic to reduce the running time. RapidNJ can be used to reconstruct large trees using a very small amount of memory by utilizing the HDD as storage." << endl << endl;
  cerr << "USAGE: rapidnj INPUT [OPTIONS]" << endl;
  cerr << "The INPUT can be a distance matrix in phylip (.phylip) format or a multiple alignment in stockholm (.sth) or phylip format (.phylip)." << endl;
  cerr << "OPTIONS:" << endl;
  cerr << "  -h, --help                display this help message and exit." << endl;
  cerr << "  -v, --verbose             turn on verbose output." << endl;
  cerr << "  -i, --input-format ARG    Specifies the type of input. pd = distance" << endl;
  cerr << "                            matrix in phylip format, sth = multiple alignment in (single line) stockholm format." << endl;
  cerr << "                            fa = multiple alignment in (single line) FASTA format." << endl;
  cerr << "  -o, --output-format ARG   Specifies the type of output. t = phylogenetic tree in newick format" << endl;
  cerr << "                            (default), m = distance matrix." << endl;
  cerr << "  -a, --evolution-model ARG Specifies which sequence evolution method to use when computing" << endl;
  cerr << "                            distance estimates from multiple alignments. jc = juke cantor," << endl;
  cerr << "                            kim = Kimura's distance (default)." << endl;
  cerr << "  -m, --memory-size         The maximum amount of memory which rapidNJ is allowed to use (in MB)." << endl;
  cerr << "                            Default is 90% of all available memory." << endl;
  cerr << "  -k, --rapidnj-mem ARG     Force RapidNJ to use a memory efficient version of rapidNJ. The 'arg'" << endl;
  cerr << "                            specifies the percentage of a sorted distance matrix which should be" << endl;
  cerr << "                            stored in memory (arg=10 means 10%)." << endl;
  cerr << "  -d, --rapidnj-disk ARG    Force RapidNJ to use HDD caching where 'arg' is the directory used to" << endl;
  cerr << "                            store cached files." << endl;
  //cerr << "  -s, --simplenj            Use a naive implementation of the NJ method." << endl;
//  cerr << "  -f, --no-rapiddist        Disable rapid computation of distance estimates and use a naive" << endl;
//  cerr << "                            algorithm for this." << endl;
  cerr << "  -c, --cores ARG           Number of cores to use for computating distance matrices from multiple" << endl;
  cerr << "                            alignments. All available cores are used by default." << endl;
  cerr << "  -b  --bootstrap ARG       Compute bootstrap values using ARG samples. The output tree will be" << endl;
  cerr << "                            annotated with the bootstrap values." << endl;
  cerr << "  -t, --alignment-type ARG  Force the input alignment to be treated as: p = protein alignment, " << endl;
  cerr << "                            d = DNA alignment." << endl;
  //cerr << "  -g  --gpu                 Use CUDA enabled GPU to compute distance estimates." << endl;
  cerr << "  -n  --no-negative-length  Adjust for negative branch lengths." << endl;
  cerr << "  -x  --output-file ARG     Output the result to this file instead of stdout." << endl;
  exit(1);
}

int main( int argc, char* argv[] ) {    

  using namespace GetOpt;
  GetOpt_pp opts(argc, argv);

  if(argc == 1 || opts >> OptionPresent('h',"help")){
    printUsage();    
  }

  opts >> OptionPresent('v', "verbose", options::verbose);
  opts >> Option('i', "input-format", options::inputFormat, "");
  opts >> Option('o', "output-format", options::outputFormat, "");
  opts >> Option('a',"evolution-model", options::distMethod, "");
  opts >> Option('m',"memory-size", options::memSize, -1);
  opts >> OptionPresent('r',"rapidnj", options::rapidNJ);
  opts >> Option('k',"rapidnj-mem", options::percentageMemoryUsage, "");
  opts >> Option('d',"rapidnj-disk", options::cacheDir, "");
  opts >> OptionPresent('s',"simplenj", options::simpleNJ);
  opts >> OptionPresent('f', "no-rapiddist", options::fastdist);
  opts >> Option('c',"cores", options::cores, -1);
  opts >> Option('b',"bootstrap", options::replicates, -1);
  opts >> Option('t',"alignment-type", options::inputtype, "");
  opts >> OptionPresent('g', "gpu", options::gpu);
  opts >> OptionPresent('n', "no-negative-length", options::negative_branches);
  opts >> Option('x', "outfile", options::outputFile, "");

  options::fastdist = !options::fastdist;
  vector<string> fileNames;
  opts >> GlobalOption(fileNames);
  if(fileNames.size() < 1){
    cerr << "ERROR: An input file must be specified!" << endl;
    exit(1);
  }
  if(fileNames.size() > 1){
    cerr << "ERROR: Only one input file can be specified!" << endl;
    exit(1);
  }
  options::fileName = fileNames[0];

  if(opts.options_remain()){
    cerr << "ERROR: One or more options were not recognised!" << endl;
    exit(0);      
  }

  streambuf* sbuf;  
  ofstream fileOutStream;

  if(options::outputFile != "") {	  
    fileOutStream.open(options::outputFile.data(), fstream::out | fstream::in | fstream::binary | fstream::trunc);
    if(!fileOutStream.is_open()){      
      cerr << "Could not open the output file \"" << options::outputFile << "\"\n";
      exit(1);
    }
    sbuf = fileOutStream.rdbuf();
  } else {
    sbuf = cout.rdbuf();
  }

  ostream out(sbuf);

  if(options::verbose){
    cerr << "RapidNJ v. " << VERSION << endl;
    if(sizeof(void*) == 4) {
      cerr << "32 bit system detected." << endl;
    } else {
      cerr << "64 bit system detected." << endl;
    }
  }

  configureNumberOfCores();

#ifdef __linux__
  gettimeofday(&startTime, NULL);
#endif

  if(options::inputFormat == ""){
    // try to determine the input format
    if (options::fileName.compare(options::fileName.length() - 4, 4, ".sth") == 0 ||
        options::fileName.compare(options::fileName.length() - 4, 4, ".STH") == 0){
      if(options::verbose){
        cerr << "Input format determined as stockholm" << endl;
      }
      options::inputFormat="sth";
    } else if (options::fileName.compare(options::fileName.length() - 3, 3, ".fa") == 0 ||
      options::fileName.compare(options::fileName.length() - 4, 4, ".faa") == 0 ||
      options::fileName.compare(options::fileName.length() - 6, 6, ".fasta") == 0 ||
      options::fileName.compare(options::fileName.length() - 3, 3, ".FA") == 0 ||
      options::fileName.compare(options::fileName.length() - 4, 4, ".FAA") == 0 ||
      options::fileName.compare(options::fileName.length() - 6, 6, ".FASTA") == 0) {
        if(options::verbose) {
          cerr << "Input format determined as FASTA" << endl;
        }
        options::inputFormat="fa";
    } else if (options::fileName.compare(options::fileName.length() - 4, 4, ".phy") == 0 || 
      options::fileName.compare(options::fileName.length() - 4, 4, ".PHY") == 0 ||
      options::fileName.compare(options::fileName.length() - 3, 3, ".ph") == 0 ||
      options::fileName.compare(options::fileName.length() - 3, 3, ".PH") == 0 ||
      options::fileName.compare(options::fileName.length() - 7, 7, ".phylip") == 0 ||
      options::fileName.compare(options::fileName.length() - 7, 7, ".PHYLIP") == 0){
      options::inputFormat = guessPhylipType(options::fileName);
      if(options::verbose){
        cerr << "Input format determined as phylip ";
        if(options::inputFormat == "pa"){
          cerr << "alignment" << endl;
        } else {
          cerr << "distance matrix" << endl;
        }
      }
    } else {
      cerr << "ERROR: could not determine input file format. Please use the '-i' option to indicate the format." << endl;
      exit(1);
    }
  }

  dataloader* dl = NULL;
  ProgressBar* pb = new ProgressBar();
    
  if(options::inputFormat == "pd"){
    //input is a distance matrix.
    if(options::replicates > -1) {
      cerr << "ERROR: Cannot perform bootstrapping with a distance matrix input." << endl;
      exit(1);
    }
    distanceMatrixInput = true;
    matrixSize = getMatrixSize();
  } else {
    dl = loadAlignment();
    matrixSize = dl->getSequenceCount();
    distanceMatrixInput = false;        
  }
  
  if(options::outputFormat == "m") {
    double memSize = getMemSize();
    double requiredMemory = 0;
    if(dl == NULL) {
      cerr << "ERROR: Both input and output format is a distance matrix." << endl;
      exit(1);
    }
    if(dl->type == DNA) {
      requiredMemory += matrixSize * dl->getSequenceLength() / 2.0;
    } else {
      requiredMemory += matrixSize * dl->getSequenceLength();
    }    
    requiredMemory += matrixSize * matrixSize * (double)sizeof(distType);    
    if(memSize < requiredMemory || options::cacheDir != "") {
      //use a disk matrix
      if(options::verbose) {
        cerr << "Using disk based matrix" << endl;
      }
      computeDistanceMatrix(true, out, true, dl);  
    } else {
      computeDistanceMatrix(false, out, true, dl);  
    }    
  } else {
    if(options::replicates > -1) {
      pb->childProgress(1.0 / (options::replicates + 1.0));
    }
    polytree* tree = computeTree(out, dl, pb);
    if(options::replicates > -1) {
      bootstrapTree(out, tree, dl, pb);
      tree->serialize_tree(out);
    } else {    
      tree->serialize_tree(out);
    }
    delete tree;
  }

  if(!distanceMatrixInput) {
    delete dl;
  }
  delete pb;

  if(options::verbose){
    cerr << endl;
    //printTime("Total running time:");
  }

  out.flush();
  if(fileOutStream.is_open()) {
    fileOutStream.close();
  }
  return 0;
}
