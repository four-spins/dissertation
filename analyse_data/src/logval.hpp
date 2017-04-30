# pragma once
# include <cmath>
# include <cassert>
# include <algorithm>
# include <sstream>
# include <iomanip>
# include <iostream>
# include <cerrno>

# ifdef COMPLEX_LOGVAL_MAKE_EIGEN_EXTENSION
# include <Eigen/Core>
# endif

namespace logval_detail {
	static double rtol = 1e-10; // default relative tolerance for comparisons
	static double atol = 1e-10; // default absolute tolerance for comparisons
	static int logval_xalloc = std::ios_base::xalloc();
	const long print_double = 0;
	const long print_native = 1;
}


std::ios_base& logval_double(std::ios_base& os)
{
  os.iword(logval_detail::logval_xalloc) = 0;
  return os;
}
 
std::ios_base& logval_native(std::ios_base& os)
{
  os.iword(logval_detail::logval_xalloc) = 1;
  return os;
}


class logval {
	// Represents a double x, which is stored logarithmically.
	// Needed for very large numbers (e.g. partition functions).
	// Provides all necessary operators to do real arithmetic
	// by only using the stored logarithms.
private:
	char sign; // -1 for negative, +1 for positive and 0 for ZERO
	double lnx; // the logarithm of x

	double rtol; // relative tolerance for comparison
	double atol; // absolute tolerance for comparison

	double log_minus(const double, const double) const; // logarithmic substraction
	double log_plus(const double, const double) const; // logarithmic addition

public:

  // converts a double to its logval expression
	logval(const double x = 0) 
    : sign {0}, lnx {x == 0.0 ? -1 : log(std::fabs(x))},
    rtol{logval_detail::rtol}, atol{logval_detail::atol}
  {
    if (x > 0) sign = +1;
    else if (x < 0) sign = -1;
    else sign = 0;
  }

	logval & operator=(const logval& rhs);
	logval & operator+=(const logval& rhs);
	logval & operator-=(const logval& rhs);
	logval & operator*=(const logval& rhs);
	logval & operator/=(const logval& rhs);

	void set_rtol(double rtol_) { rtol = rtol_; }
	void set_atol(double atol_) { atol = atol_; }

	bool operator==(const logval& rhs) const;
	bool operator!=(const logval& rhs) const;
	bool operator<(const logval& rhs) const;
	bool operator>(const logval& rhs) const;
	bool operator<=(const logval& rhs) const;
	bool operator>=(const logval& rhs) const;

	logval operator-() const;

  explicit operator double() const { return this->sign * exp(this->lnx); }
  double get_exponent() const { return this->lnx; }
  double get_sign() const { return this->sign; }

  friend inline bool isfinite(const logval&);
#ifdef COMPLEX_LOGVAL_MAKE_EIGEN_EXTENSION
	explicit operator int() const { return static_cast<int>(this->sign*exp(this->lnx)); }
#endif

	friend logval make_logval_from(int, double); // creates logval directly from sign and exponent
  std::ostream& print( std::ostream& os) const; // prints the logval as regular double
  std::ostream& print_native( std::ostream& os) const; // prints logval in 'native' format, i.e sign and exponent
	
};

std::ostream & operator<<(std::ostream &os, const logval& p) {
	if (os.iword(logval_detail::logval_xalloc) == logval_detail::print_double){
		return p.print(os);
	}
  else {
		return p.print_native(os);
  }
}

double logval::log_plus(double a, double b) const {
	// for numerical precision, it is better to have
	// the exponent in the exponential be negative
	// (harmless underflows instead of harmful overflows)
	if (b >= a) {
		return b + log1p(exp(a - b));
	}
	else {
		return a + log1p(exp(b - a));
	}
}

double logval::log_minus(double a, double b) const {
	
	// calculates the logarithm c=log(c_) of c_ = a_ - b_, inferred only from
	// the a = log(a_) and b = log(b_)
	// for numerical precision, we calculate the substraction
	// differently, depending on the magnitude of a and b
	double result = a + log1p(-exp(b - a));
	assert(a >= b);
	//if (b >= a) {
	//	// log (a-b) = log(b(a/b - 1) = log(b)
	//	return b + log(expm1(a - b));
	//}
	//else {
	return result;
	//}
}

logval& logval::operator=(const logval& rhs) {
	this->sign = rhs.sign;
	this->lnx = rhs.lnx;
	return *this;
}

logval& logval::operator+=(const logval& rhs) {
	//double bigger = std::max(this->lnx, rhs.lnx);
	//double smaller = std::min(this->lnx, rhs.lnx);
	if (0 == this->sign) {
		*this = rhs;
		return *this;
	}
	if (0 == rhs.sign) { // don't change anything if adding a zero
		return *this;
	} 

	const char S = this->sign*rhs.sign;
	// if the operands have different signs, we have
	// a subtraction
	if (S < 0) {
		if (this->sign < 0) { // first operand is negative, second positive
			// calculate (-a)+b = -a + b = b - a instead
			if (rhs.lnx > this->lnx) { 
				this->sign = +1; 
				this->lnx = log_minus(rhs.lnx, this->lnx);
			}
			else if (rhs.lnx < this->lnx) { 
				this->sign = -1;
				this->lnx = log_minus(this->lnx, rhs.lnx);
			}
			else { this->sign = 0; }
			//logval result = rhs;
			//result -= *this;
			//return result;
		}
		else { // first operand is positive, second negative
			if (this->lnx > rhs.lnx) { 
				this->sign = +1; 
				this->lnx = log_minus(this->lnx, rhs.lnx);
			}
			else if (this->lnx < rhs.lnx) { 
				this->sign = -1; 
				this->lnx = log_minus(rhs.lnx, this->lnx);
			}
			else { this->sign = 0; }
			//logval result = *this;
			//result -= rhs;
			//return result;
		}
	}
	// if both operands have equal sign, we really 
	// have an addition
	if (S > 0) {
		if (this->sign < 0) { // both operands are negative
			// calculate (-a) + (-b) = -(a+b)
			double lnx = log_plus(this->lnx, rhs.lnx);
			this->lnx = lnx;
		}
		else { // both operands are positive
			double lnx = log_plus(this->lnx, rhs.lnx);
			this->lnx = lnx;
			this->sign = +1;
		}
	}
	
	return *this;
}
logval& logval::operator-=(const logval& rhs) {
	*this += (-rhs);
	return *this;
}
logval& logval::operator*=(const logval& rhs) {
	this->lnx += rhs.lnx;
	this->sign *= rhs.sign;
	return *this;
}
logval& logval::operator/=(const logval& rhs) {
	this->lnx -= rhs.lnx;
	this->sign *= rhs.sign;
	return *this;
}

bool logval::operator==(const logval& rhs) const {
	if (this->sign != rhs.sign) { return false; }
	
	if (this->sign == 0) { return true; } // Zero is always equal
	// else: we need to compare the value of lnx...

	// sign matches, hence compare the logarithms up to 
	// some numerical deviations
	// get the smaller tolerances from both operands:
	const double atol = std::min(this->atol, rhs.atol);
	const double rtol = std::min(this->rtol, rhs.rtol);
	// compare (asymmetrically!)
	if (std::fabs(this->lnx - rhs.lnx) <= (atol + rtol * std::fabs(rhs.lnx))) {
		return true;
	}
	return false;
}

bool logval::operator!=(const logval& rhs) const {
	return !(*this == rhs);
}


logval operator+(logval lhs, const logval& rhs) {
	lhs += rhs;
	return lhs;
}
logval operator-(logval lhs, const logval& rhs) {
	lhs -= rhs;
	return lhs;
}
logval operator*(logval lhs, const logval& rhs) {
	lhs *= rhs;
	return lhs;
}
logval operator/(logval lhs, const logval& rhs) {
	lhs /= rhs;
	return lhs;
}

logval logval::operator-() const {
	logval result = *this;
	result.sign = -this->sign;
	return result;
}
bool logval::operator<(const logval& rhs) const {
	// comparisons can be made by comparing the logarithms,
	// because log(x) is monotonic
	if (this->sign == -1 && rhs.sign == -1) { return this->lnx > rhs.lnx; }
	else if (this->sign < rhs.sign) { return true; }
	else if (this->sign == +1 && rhs.sign == +1) { return this->lnx < rhs.lnx; }
	else return false;
}
bool logval::operator>(const logval& rhs) const {
	return !(*this <= rhs);
}

bool logval::operator<=(const logval& rhs) const {
	return (*this < rhs || *this == rhs);
}

bool logval::operator>=(const logval& rhs) const {
	return !(*this < rhs);
}

std::ostream& logval::print_native(std::ostream& os) const {
	if (this->sign == 0) { os << "0"; return os; }
	else if (this->sign == -1) { os << "-exp("; }
	else { os << "+exp("; }
	os << this->lnx << ")";
	return os;
}

std::ostream& logval::print(std::ostream& os) const {
	os << this->sign * exp(this->lnx);
	return os;
}

logval pow(const logval& x, unsigned p) {
	logval result = 1;
	for (unsigned i = 0; i < p; ++i) {
		result *= x;
	}
	return result;
}

logval make_logval_from(int sign, double lnx) {
	// creates logval directly from sign and exponent
	assert(-1 == sign || 1 == sign || 0 == sign);
	logval result;
	result.lnx = lnx;
	result.sign = sign;
	return result;
}

inline logval log(const logval& x) {
  errno = 0;
  if ( x.get_sign() <= 0 ){
    errno = EDOM;
    return logval {0.};
  }
  return logval(x.get_exponent());
}

inline logval sqrt(const logval& x) {
  errno = 0;
  if ( x.get_sign() < 0 ){
    errno = EDOM;
    return logval {0.};
  }
  return make_logval_from(x.get_sign(), x.get_exponent() / 2.0);
}

inline logval exp(const logval& x) {
  return make_logval_from(+1, x.get_sign()*exp(x.get_exponent()));
}

std::istream & operator>>(std::istream &is, logval& x) {
	if (is.iword(logval_detail::logval_xalloc) == logval_detail::print_double){
    // todo: exception handling
    double tmp;
    is >> tmp;
    x = logval(tmp); 
	}
  else {
    std::string val;
    is >> val;
    size_t pos_exp = val.find_first_of("exp(");
    size_t pos_exp_end = val.find_first_of(")");
    if ( std::string::npos == pos_exp ){
      x = logval(0);
    }
    else{
      double tmp = std::stod(val.substr(pos_exp+4, pos_exp_end));
      if (val.substr(0, pos_exp) == "+"){
        x = make_logval_from(+1, tmp); 
      }
      else { // negative sign
        x = make_logval_from(-1, tmp); 
      }
    }
  }
  return is;
}

inline logval abs(const logval&  x)  {
  return make_logval_from(std::abs(x.get_sign()), x.get_exponent()); 
}

inline bool isfinite(const logval&  x)  {
  return std::isfinite(x.lnx); 
}

inline logval sinh(const logval& x){
  return 0.5*(exp(x) - exp(-x));
}

inline logval cosh(const logval& x){
  return 0.5*(exp(x) + exp(-x));
}

#ifdef COMPLEX_LOGVAL_MAKE_EIGEN_EXTENSION
namespace Eigen {
	template<> struct NumTraits<logval> {
		typedef logval Real;
		typedef logval NonInteger;
		typedef logval Nested;

		enum {
			IsComplex = 0,
			IsInteger = 0,
			ReadCost = 1,
			AddCost = 10,
			MulCost = 2, //that's right, multiplying logvals is cheaper than adding them
			IsSigned = 1,
			RequireInitialization = 0
		};

		// implicit casts
		static logval epsilon() { return logval(std::numeric_limits<double>::epsilon()); }
		logval dummy_precision() { return 1e-14; }

	};
}
#endif
