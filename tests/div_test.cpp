#include <iostream>
#include "src/include/int2048.h"
using namespace sjtu;void P(const char* tag,int2048 &x){std::cout<<tag<<":"; x.print(); std::cout<<"\n";}
int main(){int2048 a; a.read("10"); int2048 b; b.read("3"); int2048 c; c.read("-10"); int2048 d; d.read("-3"); auto q1=a/b; auto q2=c/b; auto q3=a/d; auto q4=c/d; P("10/3",q1); P("-10/3",q2); P("10/-3",q3); P("-10/-3",q4);
 auto r1=a%b; auto r2=c%b; auto r3=a%d; auto r4=c%d; P("10%3", r1); P("-10%3", r2); P("10%-3", r3); P("-10%-3", r4);
 return 0;}
