#include <bits/stdc++.h>
#include "src/include/int2048.h"
using namespace sjtu;
int main(){
  int2048 a; a.read("12345678901234567890");
  int2048 b; b.read("9876543210");
  int2048 c=add(a,b);
  c.print(); std::cout<<"\n"; // expect 12345688777777780100
  int2048 d=minus(a,b);
  d.print(); std::cout<<"\n"; // expect 12345669024691345680
  int2048 e; e.read("-10"); int2048 f; f.read("3");
  add(e,f).print(); std::cout<<"\n"; // -7
  minus(e,f).print(); std::cout<<"\n"; // -13
  int2048 g; g.read("1000"); int2048 h; h.read("1000"); minus(g,h).print(); std::cout<<"\n"; // 0
  return 0;
}

