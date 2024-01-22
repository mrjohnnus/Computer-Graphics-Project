#include <cmath>

struct Vec3
{
    float x,y,z;
};

struct mat3x3
{
    float mat[3][3];
};

struct Coeff
{
    float alpha, beta, gamma;
};

Vec3 vec_sub(Vec3 v1, Vec3 v2)
{
    return {v1.x - v2.x, v1.y - v2.y, v1.z - v2.z};
}

Vec3 vec_mul(Vec3 v1, float k)
{
    return {v1.x * k, v1.y * k, v1.z * k};
}

Vec3 vec_div(Vec3 v1, float k)
{
    return {v1.x / k, v1.y / k, v1.z / k};
}

float vec_dotProduct(Vec3 v1, Vec3 v2)
{
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

Vec3 ortoghonalize(Vec3 v1, Vec3 v2)
{
    float k = vec_dotProduct(v1, v2)/vec_dotProduct(v2,v2);
    Vec3 vec_mult = vec_mul(v2, k);
    return vec_sub(v1, vec_mult);
}

Vec3 vec_crossProduct(Vec3 v1, Vec3 v2)
{
    Vec3 v;
    v.x = v1.y * v2.z - v1.z * v2.y;
    v.y = v1.z * v2.x - v1.x * v2.z;
    v.z = v1.x * v2.y - v1.y * v2.x;
    return v;
}

Vec3 vec_normalise(Vec3 v)
{
	float l = std::sqrt(vec_dotProduct(v, v));
    float x, y, z;
    x = (v.x / l) + 0.0;
    y = (v.y / l) + 0.0;
    z = (v.z / l ) + 0.0;
	return { x, y, z};
}

mat3x3 change_basis_matrix(Vec3 u, Vec3 v, Vec3 n)
{
    return {u.x, u.y, u.z,
            v.x, v.y, v.z,
            n.x, n.y, n.z};
}

Vec3 matrix_multVec(mat3x3 m, Vec3 i)
{    
    return {i.x * m.mat[0][0] + i.y * m.mat[0][1] + i.z * m.mat[0][2],
            i.x * m.mat[1][0] + i.y * m.mat[1][1] + i.z * m.mat[1][2],
            i.x * m.mat[2][0] + i.y * m.mat[2][1] + i.z * m.mat[2][2]};
}

Coeff barycentric_coeff(float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3)
{
    float t[2][2] = {{x1 - x3, x2 - x3},{y1 - y3, y2 - y3}};
    float k = 1/(t[0][0]*t[1][1] - t[0][1]*t[1][0]);
    float t_inv[2][2] = {{k*t[1][1], k*(-t[0][1])},{k*(-t[1][0]), k*t[0][0]}};
    
    float alpha, beta, gamma;
    
    alpha = t_inv[0][0]*(x0 - x3) + t_inv[0][1]*(y0 - y3);
    beta = t_inv[1][0]*(x0 - x3) + t_inv[1][1]*(y0 - y3);
    gamma = 1 - alpha - beta;

    return (Coeff) {alpha,beta,gamma};
}

Vec3 component_mul(Vec3 v1, Vec3 v2, Vec3 v3)
{
    return {v1.x*v2.x*v3.x, v1.y*v2.y*v3.y, v1.z*v2.z*v3.z};
}