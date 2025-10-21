#include <bits/stdc++.h>
#include "src/include/int2048.h"
using namespace sjtu;int main(){
  auto P=[&](const int2048 &x){x.print(); std::cout<<"\n";};
  int2048 a; a.read("5"); int2048 b; b.read("-5"); P(add(a,b)); // 0
  int2048 c; c.read("-5"); int2048 d; d.read("-5"); P(minus(c,d)); // 0
  int2048 e; e.read("0000"); P(e); // 0
  int2048 f; f.read("+10"); int2048 g; g.read("3"); P(add(f,g)); // 13
  int2048 h; h.read("-3"); P(minus(f,h)); // 13
  return 0;}

