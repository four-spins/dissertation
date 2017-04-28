/*
 * logval.h
 *
 *  Created on: Jan 12, 2011
 *      Author: gerlach
 */

#ifndef LOGVAL_H_
#define LOGVAL_H_

#include <cmath>
#include <iostream>
//calculations with logval's internally use the natural logarithm of values
//to avoid overflows
class logval {
public:
	double lnx;

	logval(); //init to 0 (lnx = -infinity)
	explicit logval(double x);
	logval(const logval& lv);

	logval& addSmaller(logval smaller); //add a smaller number
	logval& addLarger(logval larger); //add a larger number
	logval& subtractSmaller(logval smaller);
	logval& subtractLarger(logval larger);

	logval& operator *=(logval rhs);
	logval& operator /=(logval rhs);
	logval& operator +=(logval rhs); //these use an additional check as compared to addSmaller and addLarger
	logval& operator -=(logval rhs); //of course only for rhs < lhs!

	logval& operator *=(double rhs);
	logval& operator /=(double rhs);
    logval& operator ++();    
    logval  operator ++(int unused);

	bool operator==(logval rhs) const;
	bool operator!=(logval rhs) const;
	bool operator<=(logval rhs) const;
	bool operator>=(logval rhs) const;
	bool operator<(logval rhs) const;
	bool operator>(logval rhs) const;

    double toDouble() const;
    double value() const;
};

// return the actual value
inline double logval::value() const { return lnx; }

inline double logval::toDouble() const { return std::exp(lnx); }

/*
logval tologvalExp(double exponent);	//returns a logval representing exp(exponent), i.e. set lnx to exponent

logval operator+(logval a, logval b);
logval operator-(logval a, logval b);
logval operator*(logval a, logval b);
logval operator/(logval a, logval b);
logval operator*(logval a, double b);
logval operator/(logval a, double b);

logval pow(logval base, double exponent);

std::ostream& operator<<(std::ostream& stream, logval val);
*/

inline bool logval::operator==(logval rhs) const {
	return lnx == rhs.lnx;
}
inline bool logval::operator!=(logval rhs) const {
	return lnx != rhs.lnx;
}
inline bool logval::operator<=(logval rhs) const {
	return lnx <= rhs.lnx;
}
inline bool logval::operator>=(logval rhs) const {
	return lnx >= rhs.lnx;
}
inline bool logval::operator<(logval rhs) const {
	return lnx < rhs.lnx;
}
inline bool logval::operator>(logval rhs) const {
	return lnx > rhs.lnx;
}

inline logval::logval() {
	lnx = std::log(0);
}

inline logval::logval(double x) {
	lnx = std::log(x);
}

inline logval::logval(const logval& lv) {
	lnx = lv.lnx;
}

inline logval& logval::addSmaller(logval smaller) {
	lnx += log1p(std::exp(smaller.lnx - lnx));
	return *this;
}

inline logval& logval::subtractSmaller(logval smaller) {
	lnx += log1p(-std::exp(smaller.lnx - lnx));
	return *this;
}

inline logval& logval::addLarger(logval larger) {
	lnx = larger.lnx + log1p(std::exp(lnx - larger.lnx));
	return *this;
}

inline logval& logval::subtractLarger(logval larger) {
	lnx = larger.lnx + log1p(-std::exp(lnx - larger.lnx));
	return *this;
}

inline logval& logval::operator ++()
{
    logval temp(1.0);
    *this += temp;
    return *this;
}

inline logval logval::operator ++(int)
{
    logval result = *this;
    ++(*this);
    return result;
}

inline logval& logval::operator *=(logval rhs) {
	lnx += rhs.lnx;
	return *this;
}

inline logval& logval::operator /=(logval rhs) {
	lnx -= rhs.lnx;
	return *this;
}

inline logval& logval::operator +=(logval rhs) {
	if (lnx == std::log(0)) {
		lnx = rhs.lnx;
	} else if (rhs.lnx <= lnx) {
		addSmaller(rhs);
	} else {
		addLarger(rhs);
	}
	return *this;
}

inline logval& logval::operator -=(logval rhs) {
	if (rhs.lnx <= lnx) {
		subtractSmaller(rhs);
	} else {
		subtractLarger(rhs);
	}
	return *this;
}

inline logval& logval::operator *=(double rhs) {
	lnx += std::log(rhs);
	return *this;
}

inline logval& logval::operator /=(double rhs) {
	lnx -= std::log(rhs);
	return *this;
}

inline logval operator+(logval a, logval b) {
	logval t = a;
	return t += b;
}
inline logval operator-(logval a, logval b) {
	logval t = a;
	return t -= b;
}
inline logval operator*(logval a, logval b) {
	logval t = a;
	return t *= b;
}
inline logval operator/(logval a, logval b) {
	logval t = a;
	return t /= b;
}
inline logval operator*(logval a, double b) {
	logval t = a;
	return t *= b;
}
inline logval operator/(logval a, double b) {
	logval t = a;
	return t /= b;
}

inline logval pow(logval base, double exponent) {
	logval r;
	r.lnx = base.lnx * exponent;
	return r;
}

inline std::ostream& operator<<(std::ostream& stream, logval val) {
//	return stream << std::exp(val.lnx);
//	return stream << "exp(" << val.lnx << ")";
	return stream << val.lnx;
}

inline std::istream& operator>>(std::istream& stream, logval& val) {
//	return stream << std::exp(val.lnx);
//	return stream << "exp(" << val.lnx << ")";

	return stream >> val.lnx;
}

inline logval tologvalExp(double exponent) {
	logval temp;
	temp.lnx = exponent;
	return temp;
}

#endif /* LOGVAL_H_ */
