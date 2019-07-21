#include <TChain.h>
#include <vector>

constexpr int MinPointsPerEvent = 3; // reject event if N points too small

struct Point {
  static constexpr float PitchX = 0.03, PitchY = 0.05, PitchZ = 0.03;  
  //Declaration of leaves types
  Int_t           PixelNb;
  Int_t           TrackID;
  Int_t           PartPDGCode;
  Int_t           SensorNb;
  Double_t        Edep;
  Int_t           EventID;
  Double_t        PixelCentreXPosition;
  Double_t        PixelCentreYPosition;
  Double_t        PixelCentreZPosition;
  Int_t           BigSensorNb;
  void print() const {
    int ix(0),iy(0),iz(0);
    xyz2cell(ix,iy,iz);
    printf("Ev:%4d XYZ:{%+.3e,%+.3e,%+.3e}={%4d/%4d/%4d} Tr:%3d PDG:%+4d Sn:%2d Ed:%.2e BigSn:%2d Pix:%d\n",
	   EventID,PixelCentreXPosition,PixelCentreYPosition,PixelCentreZPosition,ix,iy,iz,TrackID,PartPDGCode,
	   SensorNb,Edep,BigSensorNb,PixelNb);
  }
  
  void xyz2cell(int &ix, int &iy, int &iz) const {
    ix = PixelCentreXPosition/PitchX;
    iy = PixelCentreYPosition/PitchY;
    iz = PixelCentreZPosition/PitchZ;
  }

  int compare(const Point& pntB) const {
    // ordering in X,Y,Z
    if (PixelCentreXPosition < pntB.PixelCentreXPosition) return -1;
    if (PixelCentreXPosition > pntB.PixelCentreXPosition) return  1;
    if (PixelCentreYPosition < pntB.PixelCentreYPosition) return -1;
    if (PixelCentreYPosition > pntB.PixelCentreYPosition) return  1;
    if (PixelCentreZPosition < pntB.PixelCentreZPosition) return -1;
    if (PixelCentreZPosition > pntB.PixelCentreZPosition) return  1;
    return 0;
  }

  bool isNeighbour(const Point& pntB) const {
    // check if points are neighbours (difference in pixels <=1 in all dimensions
    if (std::abs(PixelCentreXPosition - pntB.PixelCentreXPosition)>1.5*PitchX ||
	std::abs(PixelCentreYPosition - pntB.PixelCentreYPosition)>1.5*PitchY ||
	std::abs(PixelCentreZPosition - pntB.PixelCentreZPosition)>1.5*PitchZ) {
      return false;
    }
    else {
      return true;
    }
      
  }
  
  ClassDefNV(Point,1);
};


TChain inpTree;
Point point, *pointPtr = &point;
Long64_t entID = -1;
Long64_t currEventID = -1;

std::vector<Point> pointsPool;
std::vector<int> pointsPoolIdx;



struct Chain
{
  int status = 0;
  std::vector<int> pointIDs; // set of connected poits

  Chain(int pntID) {
    status = 0;
    addPoint(pntID);
  }

  void addPoint(int ip) {
    pointIDs.push_back(ip);
  }
  
  const Point& getGlobalPoint(int i) const {
    return pointsPool[ pointsPoolIdx[i] ];
  }
  
  bool isNeighbourToLast(const Point& pnt) const {
    if (!pointIDs.size()) return false; // should not happen
    const auto& last = getGlobalPoint(pointIDs.back());
    return last.isNeighbour(pnt);
  }

  void print() const {
    int np = pointIDs.size();
    printf("Stat:%d, NPnt:%3d | ",status, np);
    for (int i=0;i<np;i++) {
      printf("Pnt %3d |",i);
      getGlobalPoint(pointIDs[i]).print();
    }
  }
  
};
std::vector<Chain> chains;


int initTree(const char* nameFile = "B2_Digit_5000ev.root", const char* nameTree = "Interest")
{
  inpTree.SetName(nameTree);
  inpTree.AddFile(nameFile);
  inpTree.SetBranchAddress("PixelNb",&pointPtr->PixelNb);
  inpTree.SetBranchAddress("TrackID",&pointPtr->TrackID);
  inpTree.SetBranchAddress("PartPDGCode",&pointPtr->PartPDGCode);
  inpTree.SetBranchAddress("SensorNb",&pointPtr->SensorNb);
  inpTree.SetBranchAddress("Edep",&pointPtr->Edep);
  inpTree.SetBranchAddress("EventID",&pointPtr->EventID);
  inpTree.SetBranchAddress("PixelCentreXPosition",&pointPtr->PixelCentreXPosition);
  inpTree.SetBranchAddress("PixelCentreYPosition",&pointPtr->PixelCentreYPosition);
  inpTree.SetBranchAddress("PixelCentreZPosition",&pointPtr->PixelCentreZPosition);
  inpTree.SetBranchAddress("BigSensorNb",&pointPtr->BigSensorNb);
  inpTree.GetEntry(0);
  entID = 0;
  return inpTree.GetEntries();
}


int getNextEvent()
{
  pointsPool.clear();
  chains.clear();
  bool first = true;
  if (entID>=inpTree.GetEntries()) {
    return -1;
  }
  while(entID<inpTree.GetEntries()) {
    inpTree.GetEntry(entID++);
    if (first) {
      currEventID = point.EventID;
      first = false;
    }
    if (point.EventID!=currEventID) break;
    pointsPool.push_back(point);
  }
  printf("\nread %d points for event %lld\n", (int)pointsPool.size(), currEventID);
  return pointsPool.size();
}

void ProcessEvent()
{
  // sort points in X (beam direction) then Y then Z
  int npnt = pointsPool.size();

  if (npnt<MinPointsPerEvent) return;
  pointsPoolIdx.clear();
  pointsPoolIdx.resize(npnt);
  for (int i=npnt;i--;) pointsPoolIdx[i] = i;
  std::sort(pointsPoolIdx.begin(), pointsPoolIdx.end(), [](int a, int b) {
      int res = pointsPool[a].compare( pointsPool[b] );
      return res<0;
    });
  //
  /*
  for (int ipt=0;ipt<npnt;ipt++) {
    const auto& pnt = pointsPool[ pointsPoolIdx[ipt] ];
    printf("%3d |",ipt);
    pnt.print();
  }
  */  
  
  // create 1st seed
  chains.emplace_back(0);
  
  for (int ipt=1;ipt<npnt;ipt++) {
    const auto& pnt = pointsPool[ pointsPoolIdx[ipt] ];
    // skip repeating points
    if (pnt.compare( pointsPool[ pointsPoolIdx[ipt-1] ] )==0) {
      //    printf("skip %d\n",ipt);
      continue;
    }
    // loop over all chains and check if the point is compatible
    bool attached = false;
    for (int ich=(int)chains.size();ich--;) {
      auto& ch = chains[ich];
      if (ch.isNeighbourToLast(pnt)) {
	attached = true;
	ch.addPoint(ipt);
      }      
    }
    if (!attached) { // start new chain, since the cluster is not neighbouring with any of existing ones
      chains.emplace_back(ipt);
    }
  }
  
  // print chains
  int nch = chains.size();
  printf("\nNChains: %d\n",nch);
  for (int ich=0;ich<nch;ich++) {
    const auto& ch = chains[ich];
    printf("Ch:%2d |",ich);
    ch.print();
  }
  
}

void actar(const char* nameFile = "B2_Digit_5000ev.root", const char* nameTree = "Interest")
{
  int nent = initTree(nameFile, nameTree);
  while( getNextEvent()>=0 ) {
    ProcessEvent();
  }
}

