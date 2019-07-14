#include <TChain.h>
#include <vector>

constexpr float PitchX = 0.03, PitchY = 0.05, PitchZ = 0.03;  

struct Point {
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
  ClassDefNV(Point,1);
};

TChain inpTree;
Point point, *pointPtr = &point;
Long64_t entID = -1;
Long64_t currEventID = -1;

std::vector<Point> points;
std::vector<int> pointsIdx;

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
  points.clear();
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
    points.push_back(point);
  }
  printf("read %d points for event %lld\n", (int)points.size(), currEventID);
  return points.size();
}

void ProcessEvent()
{
  // sort points in X (beam direction) then Y then Z
  pointsIdx.clear();
  pointsIdx.resize(points.size());
  for (int i=points.size();i--;) pointsIdx[i] = i;
  std::sort(pointsIdx.begin(), pointsIdx.end(), [](int a, int b) {
      const auto& pntA = points[a];
      const auto& pntB = points[b];
      
      if (pntA.PixelCentreXPosition == pntB.PixelCentreXPosition) {
	if (pntA.PixelCentreYPosition == pntB.PixelCentreYPosition) {
	  return pntA.PixelCentreZPosition < pntB.PixelCentreZPosition;
	}
	return pntA.PixelCentreYPosition < pntB.PixelCentreYPosition;
      }
      return pntA.PixelCentreXPosition < pntB.PixelCentreXPosition;
    });
}

void actar(const char* nameFile = "B2_Digit_5000ev.root", const char* nameTree = "Interest")
{
  int nent = initTree(nameFile, nameTree);
  while( getNextEvent()>=0 ) {
    ProcessEvent();
  }
}

