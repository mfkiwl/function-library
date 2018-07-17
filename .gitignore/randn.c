//生成一个均值为mean，标准差为sig的随机数
double randn(double mean, double sig) {
    double PI = 3.14159265357;
    double a, b;
    a = ((double)rand() + 1.0) / ((double)RAND_MAX + 1.0);
    b = ((double)rand() + 1.0) / ((double)RAND_MAX + 1.0);
    return mean + sqrt(-2.0*log(a))*sin(2.0*PI*b)*sig;
}
