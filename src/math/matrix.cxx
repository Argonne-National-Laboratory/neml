#include "math/matrix.h"

#include "math/nemlmath.h"

#include <numeric>

namespace neml {

FlatVector::FlatVector(size_t n) : 
    n_(n), data_(nullptr), own_(true)
{
  data_ = new double [n];
}

FlatVector::FlatVector(size_t n, double * data) :
    n_(n), data_(data), own_(false)
{

}

FlatVector::FlatVector(const std::vector<double> input) :
    FlatVector(input.size())
{
  std::copy(input.begin(), input.end(), data_);
}     

FlatVector::FlatVector(const FlatVector & other) :
    n_(other.n()), own_(true)
{
  data_ = new double [n_];
  std::copy(other.data(), other.data() + n_, data_);
}

FlatVector::~FlatVector() 
{
  if (own_) {
    delete [] data_;
  }
  data_ = nullptr;
}

void FlatVector::copy(double * data)
{
  std::copy(data, data + n_, data_);
}

Matrix::Matrix(size_t m, size_t n) : 
    m_(m), n_(n)
{
  data_ = new double [m*n];
}

Matrix::~Matrix()
{
  delete [] data_;
  data_ = nullptr;
}

FlatVector Matrix::dot(const FlatVector & other)
{
  FlatVector res(m_);
  matvec(other, res);
  return res;
}

void Matrix::matvec(const FlatVector & other, FlatVector & res)
{
  if ((other.n() != n()) || (res.n() != m())) {
    throw std::invalid_argument("Matrix and vector sizes wrong for dot"
                                " product");
  }
  mat_vec(data_, m(), other.data_, n(), res.data_);
}

const double & Matrix::operator()(size_t i, size_t j) const
{
  return data_[CINDEX(i,j,m_)];
}

double & Matrix::operator()(size_t i, size_t j)
{
  return data_[CINDEX(i,j,m_)];
}

SquareMatrix::SquareMatrix(ParameterSet & params) :
    NEMLObject(params),
    Matrix(params.get_parameter<size_t>("m"),params.get_parameter<size_t>("m"))
{
  std::string type = params.get_parameter<std::string>("type");
  std::vector<double> data = params.get_parameter<std::vector<double>>("data");
  std::vector<size_t> blocks =
      params.get_parameter<std::vector<size_t>>("blocks"); 

  if (type == "zero") {
    std::fill(data_, data_+size(), 0);
  }
  else if (type == "identity") {
    setup_id_();
  }
  else if (type == "diagonal") {
    setup_diagonal_(data);
  }
  else if (type == "diagonal_blocks") {
    setup_diagonal_blocks_(data, blocks);
  }
  else if (type == "block") {
    setup_block_(data, blocks);
  }
  else if (type == "dense") {
    if (data.size() != size()) {
      throw std::invalid_argument("Input data does not have the right shape for"
                                  " matrix");
    }
    std::copy(data.begin(), data.end(), data_);
  }
  else {
    throw std::invalid_argument("Invalid SquareMatrix initialization type");
  }
}

std::string SquareMatrix::type()
{
  return "SquareMatrix";
}

ParameterSet SquareMatrix::parameters()
{
  ParameterSet pset(SquareMatrix::type());

  pset.add_parameter<size_t>("m");

  pset.add_optional_parameter<std::string>("type", std::string("zero"));
  pset.add_optional_parameter<std::vector<double>>("data", std::vector<double>());
  pset.add_optional_parameter<std::vector<size_t>>("blocks", std::vector<size_t>());

  return pset;
}

std::unique_ptr<NEMLObject> SquareMatrix::initialize(ParameterSet & params)
{
  return neml::make_unique<SquareMatrix>(params);
}

void SquareMatrix::setup_id_()
{
  std::fill(data_, data_+size(), 0.0);
  for (size_t i = 0; i < m(); i++) {
    data_[CINDEX(i, i, m())] = 1.0;
  }
}

void SquareMatrix::setup_diagonal_(std::vector<double> & data)
{
  if (data.size() != n()) {
    throw std::invalid_argument("For diagonal initialization data vector must"
                                " have the same size as the matrix rank");
  }
  std::fill(data_, data_+size(), 0.0);
  for (size_t i = 0; i < m(); i++) {
    data_[CINDEX(i, i, m())] = data[i];
  }
}

void SquareMatrix::setup_diagonal_blocks_(std::vector<double> & data, 
                                          std::vector<size_t> & blocks)
{
  if (blocks.size() != data.size()) {
    throw std::invalid_argument("For diagonal block initialization data vector"
                                " must have the same size as the blocks vector");
  }
  check_blocks_(blocks);

  std::fill(data_, data_+size(), 0.0);

  size_t ci = 0;
  size_t cb = 0;
  for (auto bs : blocks) {
    for (size_t i = 0; i < bs; i++) {
      data_[CINDEX(ci, ci, m())] = data[cb];
      ci++;
    }
    cb++;
  }
}

void SquareMatrix::setup_block_(std::vector<double> & data, 
                                std::vector<size_t> & blocks)
{
  if (data.size() != (blocks.size() * blocks.size())) {
    throw std::invalid_argument("For block initialization the data vector"
                                " must have a size equal to the number of blocks"
                                " squared");
  }
  check_blocks_(blocks);

  auto offsets = offsets_(blocks);

  for (size_t bi = 0; bi < blocks.size(); bi++) {
    for (size_t bj = 0; bj < blocks.size(); bj++) {
      size_t co = CINDEX(bi, bj, blocks.size());
      for (size_t i = offsets[bi]; i < offsets[bi+1]; i++) {
        for (size_t j = offsets[bj]; j < offsets[bj+1]; j++) {
          data_[CINDEX(i,j,m())] = data[co];
        }
      }
    }
  }
}

void SquareMatrix::check_blocks_(std::vector<size_t> & blocks)
{
  size_t sum = 0;
  for (auto bs : blocks) {
    sum += bs;
  }
  if (sum != n()) {
    throw std::invalid_argument("The blocks vector must have entries that"
                                " sum to the matrix rank");
  }
}

std::vector<size_t> SquareMatrix::offsets_(std::vector<size_t> & blocks)
{
  std::vector<size_t> result(blocks.size() + 1);
  result[0] = 0;
  std::partial_sum(blocks.begin(), blocks.end(), result.begin() + 1);
  return result;
}

} // namespace neml
