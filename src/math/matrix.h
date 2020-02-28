#pragma once

#include "../windows.h"

#include "../objects.h"

#include <algorithm>
#include <stdexcept>

namespace neml {

/// Dense vector class
class NEML_EXPORT FlatVector {
 public:
  FlatVector(size_t n);
  FlatVector(size_t n, double * data);
  FlatVector(const std::vector<double> input);
  virtual ~FlatVector();

  size_t n() const;
  bool owns_data() const;

  void copy(double * data);

  double * data();

 protected:
  size_t n_;
  double * data_;
  bool own_;

  friend class Matrix;
};

/// Dense matrix class
class NEML_EXPORT Matrix {
 public:
  Matrix(size_t m, size_t n);
  virtual ~Matrix();

  size_t m() const;
  size_t n() const;
  size_t size() const;

  FlatVector dot(const FlatVector & other);
  void matvec(const FlatVector & other, FlatVector & res);

  double * data();

 protected:
  size_t m_, n_;
  double * data_;
};

/// Specialized square dense matrix
//    Ways I may want to initialize:
//      zero
//      identity
//      some generic constant diagonal
//      diagonal blocks
//      generic blocks
//      fully-dense
class NEML_EXPORT SquareMatrix: public Matrix, public NEMLObject {
 public:
  SquareMatrix(size_t m, std::string type, std::vector<double> data,
               std::vector<size_t> blocks);

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

 private:
  /// Initialize an identity matrix
  void setup_id_();
  /// Initialize a general diagonal matrix
  void setup_diagonal_(std::vector<double> & data);
  /// Initialize diagonal matrix with different diagonal entries
  void setup_diagonal_blocks_(std::vector<double> & data, 
                              std::vector<size_t> & blocks);
  /// Initialize by blocks
  void setup_block_(std::vector<double> & data,
                    std::vector<size_t> & blocks);
  /// Ensure that the block definitions have the right dimensions
  void check_blocks_(std::vector<size_t> & blocks);
  /// Basically cumsum, generate offsets from sizes
  std::vector<size_t> offsets_(std::vector<size_t> & blocks);
};

static Register<SquareMatrix> regSquareMatrix;

} // namespace neml
