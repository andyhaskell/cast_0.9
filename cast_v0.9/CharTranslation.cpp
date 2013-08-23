#include <cctype>
#include <cassert>

const char half_byte_to_char[] = "01234567ACGTN-is";

unsigned char char_to_half_byte(char c) {
  if (isdigit(c)) return c-'0';
  switch(c) {
  case 'A': return 8;
  case 'C': return 9;
  case 'G': return 10;
  case 'T': return 11;
    //case 'N': return 12;
  case '-': return 13;
  case 'i': return 14;
  case 's': return 15;
  default: return 12; // convert other characters to 'N'
  }
}

