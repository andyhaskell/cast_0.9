#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

#include "CharTranslation.cpp"
#include "IOUtils.cpp"

using namespace std;

class HalfByteWriter {
  char leftover; // if odd number of chars
  fstream fout;
  int tot_script_chars; // number of script chars so far

public:
  HalfByteWriter(const char *out_file) { // opens file
    mustOpenStream(fout, out_file, ios::out|ios::binary);
    leftover = 0;
    tot_script_chars = 0;
  }
  
  int get_tot_script_chars(void) {
    return tot_script_chars;
  }

  // note this has to be able to convert from chars to 0-15
  void write(const string &edit_script) {
    tot_script_chars += edit_script.length();
    int len = (leftover != 0) + edit_script.length();
    char buf[len/2];
    int buf_pos = 0, script_pos = 0;
    if (leftover) // use leftover half-byte
      buf[buf_pos++] = (char_to_half_byte(leftover)<<4) |
	char_to_half_byte(edit_script[script_pos++]);
    for (; buf_pos < len/2; buf_pos++) {
      buf[buf_pos] = (char_to_half_byte(edit_script[script_pos])<<4) |
	char_to_half_byte(edit_script[script_pos+1]);
      script_pos += 2;
    }
    if (len%2) // odd number of characters; save leftover
      leftover = edit_script[script_pos];
    else
      leftover = 0;
    fout.write(buf, len/2);
  }

  // adds 1-2 '0' characters and flushes last 1/2-byte if necessary; closes file		
  void flush_and_close(void) {
    write("00");
    fout.close();
  }
};

