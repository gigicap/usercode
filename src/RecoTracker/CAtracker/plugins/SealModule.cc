#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "RecoTracker/TkTrackingRegions/interface/OrderedHitsGeneratorFactory.h"
#include "RecoTracker/TkTrackingRegions/interface/OrderedHitsGenerator.h"
#include "CAHitsGenerator.h"
#include "CAHitsGeneratorForDebuggers.h"


DEFINE_EDM_PLUGIN(OrderedHitsGeneratorFactory, CAHitsGenerator, "CAHitsGenerator");
DEFINE_EDM_PLUGIN(OrderedHitsGeneratorFactory, CAHitsGeneratorForDebuggers, "CAHitsGeneratorForDebuggers");
