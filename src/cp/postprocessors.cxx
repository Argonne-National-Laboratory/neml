#include "cp/postprocessors.h"

namespace neml {

CrystalPostprocessor::CrystalPostprocessor(ParameterSet & params) :
    NEMLObject(params)
{

}

PTRTwinReorientation::PTRTwinReorientation(ParameterSet & params) :
    CrystalPostprocessor(params),
    threshold_(params.get_object_parameter<Interpolate>("threshold")),
    prefix_(params.get_parameter<std::string>("prefix"))
{

}

std::string PTRTwinReorientation::type()
{
  return "PTRTwinReorientation";
}

ParameterSet PTRTwinReorientation::parameters()
{
  ParameterSet pset(PTRTwinReorientation::type());

  pset.add_parameter<NEMLObject>("threshold");
  pset.add_optional_parameter<std::string>("prefix", std::string("slip"));
  
  return pset;
}

std::unique_ptr<NEMLObject> PTRTwinReorientation::initialize(
    ParameterSet & params)
{
  return neml::make_unique<PTRTwinReorientation>(params);
}

void PTRTwinReorientation::populate_hist(const Lattice & L, 
                                            History & history) const
{
  size_t j = 0;
  for (size_t g = 0; g < L.ngroup(); g++) {
    for (size_t i = 0; i < L.nslip(g); i++) {
      if (L.slip_type(g, i) == Lattice::SlipType::Twin)
        history.add<double>("twin_fraction"+std::to_string(j));
      j++;
    }
  }
  history.add<double>("twinned");
}

void PTRTwinReorientation::init_hist(const Lattice & L,
                                        History & history) const
{
  size_t j = 0;
  for (size_t g = 0; g < L.ngroup(); g++) {
    for (size_t i = 0; i < L.nslip(g); i++) {
      if (L.slip_type(g, i) == Lattice::SlipType::Twin)
        history.get<double>("twin_fraction"+std::to_string(j)) = 0.0;
      j++;
    }
  }
  history.get<double>("twinned") = 0.0;
}

void PTRTwinReorientation::act(SingleCrystalModel & model, 
                               const Lattice & L, const double & T, 
                               const Symmetric & D,
                               const Skew & W, History & state,
                               const History & prev_state)
{
  // If "not twinned"
  if (prev_state.get<double>("twinned") < 0.5) {
    state.get<double>("twinned") = prev_state.get<double>("twinned");
    size_t j = 0;
    for (size_t g = 0; g < L.ngroup(); g++) {
      for (size_t i = 0; i < L.nslip(g); i++) {
        if (L.slip_type(g, i) == Lattice::SlipType::Twin) {
          // Get the twin fraction for this system
          double slip = state.get<double>(prefix_+std::to_string(j));
          state.get<double>("twin_fraction"+std::to_string(j)) = slip / 
              L.characteristic_shear(g,i);
          // If over the threshold do stuff
          if ((state.get<double>("twin_fraction"+std::to_string(j)) >
              threshold_->value(T)) && (state.get<double>("twinned") < 0.5)) {
            Orientation twinned = L.reorientation(g,i) * 
                model.get_active_orientation(state);
            model.set_active_orientation(state, twinned);
            state.get<double>("twinned") = 1.0;
          }
        }
        j++;
      }
    }
    // Remove accumulated twin slip
    if (state.get<double>("twinned") > 0.5) {
      j = 0;
      for (size_t g = 0; g < L.ngroup(); g++) {
        for (size_t i = 0; i < L.nslip(g); i++) {
          if (L.slip_type(g,i) == Lattice::SlipType::Twin)
            state.get<double>(prefix_+std::to_string(j)) = 0;
          j++;
        }
      }
      // Allow to "retwin"
      state.get<double>("twinned") = 0.0;
    }
  }
  else {
    state.get<double>("twinned") = prev_state.get<double>("twinned");
    size_t j = 0;
    for (size_t g = 0; g < L.ngroup(); g++) {
      for (size_t i = 0; i < L.nslip(g); i++) {
        if (L.slip_type(g, i) == Lattice::SlipType::Twin) {
          state.get<double>("twin_fraction"+std::to_string(j)) = 
              prev_state.get<double>("twin_fraction"+std::to_string(j));
        }
        j++;
      }
    }
  }
}

} // namespace neml
