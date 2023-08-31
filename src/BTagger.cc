#include <BTagger.h>

#include <cmath>

#include <Options.h>       
#include <PhysicsObjects.h>


BTagger::BTagger(Options const &options)                       
    : bTagCutLoose_{Options::NodeAs<double>(
        options.GetConfig(), {"b_tagger", "tag_threshold_loose"})},
      bTagCutMedium_{Options::NodeAs<double>(
        options.GetConfig(), {"b_tagger", "tag_threshold_medium"})},
      bTagCutTight_{Options::NodeAs<double>(
        options.GetConfig(), {"b_tagger", "tag_threshold_tight"})},
      ptCut_{Options::NodeAs<double>(
        options.GetConfig(), {"b_tagger", "min_pt"})},
      etaCut_{Options::NodeAs<double>(
        options.GetConfig(), {"b_tagger", "max_abs_eta"})}{
}
                                                               
                                                               
bool BTagger::operator()(Jet const &jet) const {               
  // return Loose(jet);
  return Medium(jet);
}                                                              


bool BTagger::Loose(Jet const &jet) const {               
  if (IsTaggable(jet)) {                                       
    return (jet.bTag >= bTagCutLoose_) ? true : false;              
  }                                                            
    return false;                                                
}                                                              


bool BTagger::Medium(Jet const &jet) const {               
  if (IsTaggable(jet)) {                                       
    return (jet.bTag >= bTagCutMedium_) ? true : false;              
  }                                                            
    return false;                                                
}                                                              


bool BTagger::Tight(Jet const &jet) const {               
  if (IsTaggable(jet)) {                                       
    return (jet.bTag >= bTagCutTight_) ? true : false;              
  }                                                            
    return false;                                                
}                                                              
                                                               
                                                               
bool BTagger::IsTaggable(Jet const &jet) const {               
  if (std::abs(jet.p4.Eta()) < etaCut_ && jet.p4.Pt() > ptCut_)
    return true;                                               
  else                                                         
    return false;                                              
}                                                              

