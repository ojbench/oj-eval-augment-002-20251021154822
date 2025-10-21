#include "int2048.h"

namespace sjtu {
static const int BASE=10000;static const int BD=4;using std::vector;using std::string;using std::istream;using std::ostream;template<typename T> using cplx=std::complex<T>;static void fft(std::vector<cplx<long double>>& a,bool inv){int n=a.size();vector<int> rev(n);for(int i=1,j=0;i<n;i++){int bit=n>>1;for(;j&bit;bit>>=1)j^=bit;j^=bit;rev[i]=j;}for(int i=0;i<n;i++)if(i<rev[i])std::swap(a[i],a[rev[i]]);const long double PI=3.14159265358979323846264338327950288L;for(int len=1;len<n;len<<=1){long double ang=(inv?-1:1)*PI/len;cplx<long double> wlen(cosl(ang),sinl(ang));for(int i=0;i<n;i+=2*len){cplx<long double> w(1,0);for(int j=0;j<len;++j){auto u=a[i+j],v=a[i+j+len]*w;a[i+j]=u+v;a[i+j+len]=u-v;w*=wlen;}}}if(inv){for(auto &x:a)x/=n;}}
static void conv(const vector<int>&A,const vector<int>&B,vector<long long>&res){int n=1;while(n<(int)A.size()+(int)B.size())n<<=1;vector<cplx<long double>> fa(n),fb(n);for(size_t i=0;i<A.size();++i)fa[i]=A[i];for(size_t i=0;i<B.size();++i)fb[i]=B[i];fft(fa,false);fft(fb,false);for(int i=0;i<n;i++)fa[i]*=fb[i];fft(fa,true);res.assign(n,0);for(int i=0;i<n;i++)res[i]=(long long)llroundl(fa[i].real());}

void int2048::trim(){while(!d.empty()&&d.back()==0)d.pop_back();if(d.empty())s=0;}
int int2048::absCmp(const int2048&a,const int2048&b){if(a.d.size()!=b.d.size())return a.d.size()<b.d.size()?-1:1;for(int i=(int)a.d.size()-1;i>=0;--i)if(a.d[i]!=b.d[i])return a.d[i]<b.d[i]?-1:1;return 0;}
int2048 int2048::addAbs(const int2048&a,const int2048&b){int2048 r;r.s=1;int n=std::max(a.d.size(),b.d.size());r.d.resize(n,0);int carry=0;for(int i=0;i<n;i++){long long cur=carry+(i<(int)a.d.size()?a.d[i]:0)+(i<(int)b.d.size()?b.d[i]:0);r.d[i]=(int)(cur%BASE);carry=(int)(cur/BASE);}if(carry)r.d.push_back(carry);return r;}
int2048 int2048::subAbs(const int2048&a,const int2048&b){int2048 r;r.s=1;r.d.resize(a.d.size());int carry=0;for(size_t i=0;i<a.d.size();++i){long long cur=(long long)a.d[i]-carry-(i<b.d.size()?b.d[i]:0);if(cur<0){cur+=BASE;carry=1;}else carry=0;r.d[i]=(int)cur;}r.trim();return r;}
int2048 int2048::mulAbs(const int2048&a,const int2048&b){int2048 r;if(a.s==0||b.s==0){r.s=0;return r;}if(a.d.size()<64||b.d.size()<64){r.s=1;r.d.assign(a.d.size()+b.d.size(),0);for(size_t i=0;i<a.d.size();++i){long long carry=0;for(size_t j=0;j<b.d.size()||carry;++j){long long cur=r.d[i+j]+carry+(long long)a.d[i]*(j<b.d.size()?b.d[j]:0);r.d[i+j]=(int)(cur%BASE);carry=cur/BASE;}}r.trim();return r;}vector<long long> tmp;conv(a.d,b.d,tmp);r.s=1;r.d.clear();r.d.reserve(tmp.size());long long carry=0;for(size_t i=0;i<tmp.size();++i){long long cur=tmp[i]+carry;int digit=(int)(cur%BASE);if(digit<0){digit+=BASE;cur-=BASE;}r.d.push_back(digit);carry=cur/BASE;}while(carry){r.d.push_back((int)(carry%BASE));carry/=BASE;}r.trim();return r;}
static int cmpAbsVec(const vector<int>&a,const vector<int>&b){if(a.size()!=b.size())return a.size()<b.size()?-1:1;for(int i=(int)a.size()-1;i>=0;--i)if(a[i]!=b[i])return a[i]<b[i]?-1:1;return 0;}
static vector<int> mulVecInt(const vector<int>&a,int m){vector<int> r;a.size();if(m==0||a.empty())return r;long long carry=0; r.resize(a.size());for(size_t i=0;i<a.size();++i){long long cur=carry+(long long)a[i]*m;r[i]=(int)(cur%BASE);carry=cur/BASE;}while(carry){r.push_back((int)(carry%BASE));carry/=BASE;}return r;}
static void subVecInPlace(vector<int>&a,const vector<int>&b){int carry=0;for(size_t i=0;i<a.size();++i){long long cur=(long long)a[i]-(i<b.size()?b[i]:0)-carry;if(cur<0){cur+=BASE;carry=1;}else carry=0;a[i]=(int)cur;}while(!a.empty()&&a.back()==0)a.pop_back();}
void int2048::divAbs(const int2048&x,const int2048&y,int2048&q,int2048&r){q.d.clear();q.s=r.s=0;r.d.clear();if(y.s==0) return; if(int2048::absCmp(x,y)<0){r=x;q.s=0;r.s=x.s;return;}vector<int> a=x.d,b=y.d; r.d.clear(); r.s=0; q.d.assign(a.size(),0); for(int i=(int)a.size()-1;i>=0;--i){ if(!r.d.empty()||a[i]){ r.d.insert(r.d.begin(),a[i]); } else { r.d.push_back(a[i]); }
 while(!r.d.empty()&&r.d.back()==0)r.d.pop_back(); int low=0,high=BASE-1,best=0; while(low<=high){int mid=(low+high)>>1; vector<int> prod=mulVecInt(b,mid); int cmp=cmpAbsVec(prod,r.d); if(cmp<=0){best=mid;low=mid+1;} else high=mid-1;} q.d[i]=best; vector<int> prod=mulVecInt(b,best); subVecInPlace(r.d,prod); }
 while(!q.d.empty()&&q.d.back()==0)q.d.pop_back(); if(q.d.empty())q.s=0; else q.s=1; if(r.d.empty())r.s=0; else r.s=1;}
int2048 int2048::fromInt(long long v){int2048 r;r.s=0;if(v==0)return r;if(v<0){r.s=-1;v=-v;}else r.s=1;while(v){r.d.push_back((int)(v%BASE));v/=BASE;}return r;}

int2048::int2048(){s=0;}
int2048::int2048(long long v){*this=fromInt(v);}
int2048::int2048(const string &str){read(str);}
int2048::int2048(const int2048 &o){d=o.d;s=o.s;}

void int2048::read(const string &str){d.clear();s=0;int n=str.size();int i=0;while(i<n&&isspace((unsigned char)str[i]))++i;int r=n-1;while(r>=i&&isspace((unsigned char)str[r]))--r;int neg=0;if(i<=r&&(str[i]=='+'||str[i]=='-')){neg=(str[i]=='-');++i;}while(i<=r&&str[i]=='0')++i;vector<int> tmp;for(int j=r;j>=i;j-=BD){int x=0;int l=std::max(i,j-BD+1);for(int k=l;k<=j;++k)x=x*10+(str[k]-'0');tmp.push_back(x);}d=tmp;trim();if(!d.empty())s=neg?-1:1;}
void int2048::print(){if(s==0){std::cout<<0;return;}if(s<0)std::cout<<'-';int n=d.size();std::cout<<d.back();for(int i=n-2;i>=0;--i){char buf[8];std::snprintf(buf,sizeof(buf),"%0*d",BD,d[i]);std::cout<<buf;}}

int2048 &int2048::add(const int2048 &o){if(o.s==0)return *this;if(s==0){*this=o;return *this;}if(s==o.s){*this=addAbs(*this,o);s=o.s;return *this;}int c=absCmp(*this,o);if(c==0){d.clear();s=0;return *this;}if(c>0){int ss=s; *this=subAbs(*this,o); s=ss;}else{*this=subAbs(o,*this);s=o.s;}return *this;}
int2048 add(int2048 a,const int2048 &b){return a.add(b);}
int2048 &int2048::minus(const int2048 &o){int2048 t=o; if(t.s) t.s=-t.s; return add(t);}
int2048 minus(int2048 a,const int2048 &b){return a.minus(b);}

int2048 int2048::operator+() const{return *this;}
int2048 int2048::operator-() const{int2048 r=*this;if(r.s)r.s=-r.s;return r;}
int2048 &int2048::operator=(const int2048 &o){if(this==&o)return *this;d=o.d;s=o.s;return *this;}

int2048 &int2048::operator+=(const int2048 &o){return add(o);}
int2048 operator+(int2048 a,const int2048 &b){return a+=b;}
int2048 &int2048::operator-=(const int2048 &o){return (*this)+=(-o);}
int2048 operator-(int2048 a,const int2048 &b){return a-=b;}

int2048 &int2048::operator*=(const int2048 &o){if(s==0||o.s==0){d.clear();s=0;return *this;}int signRes=s*o.s;int2048 A=*this,B=o;A.s=B.s=1;*this=mulAbs(A,B);s=signRes;trim();return *this;}
int2048 operator*(int2048 a,const int2048 &b){return a*=b;}

int2048 &int2048::operator/=(const int2048 &o){int sa=s,sb=o.s;if(sb==0){d.clear();s=0;return *this;}if(sa==0){return *this;}int2048 A=*this,B=o;A.s=B.s=1;int2048 q,r;divAbs(A,B,q,r);if(sa==sb){*this=q; if(!d.empty())s=1; else s=0;}
else{ if(r.s==0){*this=q; s = q.d.empty()?0:-1;} else { // floor adjust: q = -(q0+1)
 if(q.s==0){q=fromInt(1);} else {int2048 one=fromInt(1); q=addAbs(q,one);} *this=q; s=-1; }
}
return *this;}
int2048 operator/(int2048 a,const int2048 &b){return a/=b;}

int2048 &int2048::operator%=(const int2048 &o){int sa=s,sb=o.s;int2048 A=*this,B=o;A.s=B.s=1;int2048 q,r;divAbs(A,B,q,r);if(sa==sb){*this=r; s=r.s;}
else{ if(r.s==0){d.clear();s=0;} else { // r = |b| - r
 int2048 bb=subAbs(B,r); *this=bb; s=sb>0?1:-1; }
} return *this;}
int2048 operator%(int2048 a,const int2048 &b){return a%=b;}

istream &operator>>(istream &is,int2048 &x){string s;is>>s;x.read(s);return is;}
ostream &operator<<(ostream &os,const int2048 &x){if(x.s==0){os<<0;return os;}if(x.s<0)os<<'-';os<<x.d.back();for(int i=(int)x.d.size()-2;i>=0;--i){char buf[8];std::snprintf(buf,sizeof(buf),"%0*d",BD,x.d[i]);os<<buf;}return os;}

bool operator==(const int2048 &a,const int2048 &b){return a.s==b.s&&a.d==b.d;}
bool operator!=(const int2048 &a,const int2048 &b){return !(a==b);}
bool operator<(const int2048 &a,const int2048 &b){if(a.s!=b.s)return a.s<b.s; if(a.s==0)return false; int c=int2048::absCmp(a,b); return a.s>0?c<0:c>0;}
bool operator>(const int2048 &a,const int2048 &b){return b<a;}
bool operator<=(const int2048 &a,const int2048 &b){return !(b<a);}
bool operator>=(const int2048 &a,const int2048 &b){return !(a<b);}

} // namespace sjtu
