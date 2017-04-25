#ifndef CREEP_H
#define CREEP_H

#include "interpolate.h"

namespace neml {

/// Scalar creep functions in terms of effective stress and strain
class ScalarCreepRule {
  public:
   virtual int g(double seq, double eeq, double t, double T, double & g) const = 0;
   virtual int dg_ds(double seq, double eeq, double t, double T, double & dg) const = 0;
   virtual int dg_de(double seq, double eeq, double t, double T, double & dg) const = 0;
   virtual int dg_dt(double seq, double eeq, double t, double T, double & dg) const;
   virtual int dg_dT(double seq, double eeq, double t, double T, double & dg) const;
};

/// Simple power law creep
class PowerLawCreep: public ScalarCreepRule {
  public:
   PowerLawCreep(double A, double n);
   PowerLawCreep(std::shared_ptr<Interpolate> A, std::shared_ptr<Interpolate> n);

   virtual int g(double seq, double eeq, double t, double T, double & g) const;
   virtual int dg_ds(double seq, double eeq, double t, double T, double & dg) const;
   virtual int dg_de(double seq, double eeq, double t, double T, double & dg) const;

   double A(double T) const;
   double n(double T) const;

  private:
   const std::shared_ptr<const Interpolate> A_, n_;
};

/// Classical Norton-Bailey creep
class NortonBaileyCreep: public ScalarCreepRule {
 public:
  NortonBaileyCreep(double A, double m, double n);
  NortonBaileyCreep(std::shared_ptr<Interpolate> A, std::shared_ptr<Interpolate> m,
                    std::shared_ptr<Interpolate> n);

   virtual int g(double seq, double eeq, double t, double T, double & g) const;
   virtual int dg_ds(double seq, double eeq, double t, double T, double & dg) const;
   virtual int dg_de(double seq, double eeq, double t, double T, double & dg) const;

   double A(double T) const;
   double m(double T) const;
   double n(double T) const;

 private:
   const std::shared_ptr<const Interpolate> A_, m_, n_;

};

} // namespace neml

#endif // CREEP_H
