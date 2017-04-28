/**
* @brief class wrapper for random numbers using Marsaglia RNG
*
* @author Marco Mueller
* @date 26.08.2010
* @note guest student programme on parallel computing at JSC
*
*/
#ifndef __RNGCLASS_MAR_HPP__
#define __RNGCLASS_MAR_HPP__

# include <ctime>
# include <sstream>
# include <string>
# include <limits>
# include <iomanip>
# include <stdexcept>
# include "dlib/serialize.h"

class rngclass {

private:
  // random number generator initialization
  unsigned seed;
  double MARSarray[98], MARSc, MARScd, MARScm;
  int MARSi, MARSj;

public:

  rngclass() : seed(static_cast<int>(time(NULL))) {
    // calculate some seed out of the system's time,
    // if there is no seed the environment variables given
    // seed the generator
    init_RAN(seed, 23, 45, 56);
  }

  rngclass(int seed_) : seed(seed_) {
    // seed the generator
    init_RAN(seed, 23, 45, 56);
  }

  void set_seed(int seed_) {
    // calculate some seed out of the system's time,
    // seed the generator
    seed = seed_;
    init_RAN(seed, 23, 45, 56);
  }

  ~rngclass() {}

  /// draws a random number out of (0,1)
  double rand() { return RAN01(); }

  const char *get_type() { return "Marsaglia RNG"; }

  unsigned get_seed() { return seed; }

  void init_RAN(int nA1, int nA2, int nA3, int nB1) {
    /*
     initializes the global data of the
     MARSAGLIA pseudo random number generator
    */

    int mA1, mA2, mA3, mANEW, mB1, mHELP;
    int i1, i2;
    double varS, varT;

    mA1 = nA1;
    mA2 = nA2;
    mA3 = nA3;
    mB1 = nB1;
    MARSi = 97;
    MARSj = 33;

    for (i1 = 1; i1 < 98; i1++) {
      varS = 0.0;
      varT = 0.5;
      for (i2 = 1; i2 < 25; i2++) {
        mANEW = (((mA1 * mA2) % 179) * mA3) % 179;
        mA1 = mA2;
        mA2 = mA3;
        mA3 = mANEW;
        mB1 = (53 * mB1 + 1) % 169;
        mHELP = (mB1 * mANEW) % 64;
        if (mHELP > 31)
          varS += varT;
        varT *= 0.5;
      }
      MARSarray[i1] = varS;
    }

    MARSc  = 362436.0 / 16777216.0;
    MARScd = 7654321.0 / 16777216.0;
    MARScm = 16777213.0 / 16777216.0;

    return;
  }

  double RAN01(void) {
    /*
     generates a pseudo random number 0 .. +1
     following the proposal of MARSAGLIA
    */

    double ranMARS;

    ranMARS = MARSarray[MARSi] - MARSarray[MARSj];
    if (ranMARS < 0.0)
      ranMARS += 1.0;

    MARSarray[MARSi] = ranMARS;

    MARSi--;
    if (MARSi < 1)
      MARSi = 97;

    MARSj--;
    if (MARSj < 1)
      MARSj = 97;

    MARSc -= MARScd;
    if (MARSc < 0.0)
      MARSc += MARScm;

    ranMARS -= MARSc;
    if (ranMARS < 0.0)
      ranMARS += 1.0;

    return ranMARS;
  }

  void from_ascii_string_to_state(std::istream &statestr){
    //std::stringstream state(instr);
//    statestr << std::setprecision( std::numeric_limits<double>::max() );
    statestr >> seed;
    for ( size_t i = 0; i < 98; i++ ){
      statestr >> MARSarray[i];
    }
    statestr >> MARSc >> MARScd >> MARScm;
    statestr >> MARSi >> MARSj;
  }

  void save(const std::string &fileName);
  void restore(const std::string &fileName);

  friend std::ostream& operator<< (std::ostream& stream, const rngclass& rng);
  friend std::istream& operator>> (std::istream& stream, rngclass& rng);
};

std::ostream& operator<< (std::ostream& stream, const rngclass& rng){
  // make sure ostream was opened with std::ios::binary!
  dlib::serialize(rng.MARSi, stream);
  dlib::serialize(rng.MARSj, stream);
  dlib::serialize(rng.MARSc, stream);
  dlib::serialize(rng.MARScd, stream);
  dlib::serialize(rng.MARScm, stream);
  dlib::serialize(rng.MARSarray, stream);
  return stream;
}
std::istream& operator>> (std::istream& stream, rngclass& rng){
  // make sure istream was opened with std::ios::binary!
  dlib::deserialize(rng.MARSi, stream);
  dlib::deserialize(rng.MARSj, stream);
  dlib::deserialize(rng.MARSc, stream);
  dlib::deserialize(rng.MARScd, stream);
  dlib::deserialize(rng.MARScm, stream);
  dlib::deserialize(rng.MARSarray, stream);
  return stream;
}
#endif //__RNGCLASS_MAR_HPP__
