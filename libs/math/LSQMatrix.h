#ifndef mifit_math_LSQMatrix_h
#define mifit_math_LSQMatrix_h

class LSQMatrix
{

    float mat[3][3];
    float vec[3];
    bool init;

public:

    LSQMatrix();
    bool IsOK();
    void Void();
    bool Load(const char *pathname);
    bool Save(const char *pathname);
    float Xvalue(float, float, float);
    float Yvalue(float, float, float);
    float Zvalue(float, float, float);
    void SetMatrix(float mat[3][3], float v[3]);
    void SetMatrix(double mat[3][3], double v[3]);
    void GetMatrix(double mat[3][3], double v[3]);
};

inline LSQMatrix::LSQMatrix()
{
    init = 0;
}

inline bool LSQMatrix::IsOK()
{
    return init;
}

inline void LSQMatrix::Void()
{
    init = false;
}

#endif // ifndef mifit_math_LSQMatrix_h
