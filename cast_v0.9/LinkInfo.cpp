struct LinkInfo {
  int start0;       // start index on original string
  int start1;       // start index on string in uniques db
  int start_script; // start index of edit script (in scripts file)
  short dist0;      // end-start on original string
  short dist1;      // end-start on string in uniques db

  LinkInfo(int _start0=0, int _start1=0, int _start_script=0, short _dist0=0, short _dist1=0) :
    start0(_start0), start1(_start1), start_script(_start_script), dist0(_dist0), dist1(_dist1) {};
};
