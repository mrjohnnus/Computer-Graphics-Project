#include "screen.h"
#include "vecs.h"
#include <fstream>
#include <cmath>
#include <sstream>
#include <iostream>
#include <limits>
#include <unordered_map>

struct Triangle
{
        std::vector<Vec3> points;
        std::vector<Vec3> points_screen;
        Vec3 normal;
};

struct Mesh
{
        std::vector<Triangle> triangles;
        
};

struct Camera
{
        Vec3 C, N, V, U;
        float d, hx, hy;
};

struct Lighting
{
    Vec3 Iamb, Il, Pl, Kd, Od;
    float Ka, Ks, eta;
};

float zbuffer[WIDTH][HEIGHT];

bool compare(const std::pair<int, int>& a, const std::pair<int, int>& b) {
    return a.first > b.first;
}

void draw_Line(Screen& screen, int x1, int y1, int x2, int y2) 
{
        Vec3 colorr = {0,255,0};
        double x = x2 - x1; 
	double y = y2 - y1; 
	double length = sqrt( x*x + y*y ); 
	double addx = x / length; 
	double addy = y / length; 
	x = x1; 
	y = y1; 
	
	for (int i = 0; i < length; i += 1) { 
		screen.pixel(x, y,colorr.x, colorr.y, colorr.z); 
		x += addx; 
		y += addy; 
		
	} 
}

void read_file(std::vector<float>& numbers)
{
        std::ifstream file("object.byu");
        float number;
        while (file >> number)
                numbers.push_back(number);

        file.close();
}

void read_camera_file(Camera& camera)
{
        std::ifstream file("camera.txt");
        std::string line;
        while (std::getline(file, line)) 
        {
                std::istringstream iss(line.substr(line.find("=") + 1));
                std::vector <float> aux_vec;
                float aux_num;
            
                while(iss >> aux_num)
                {
                        aux_vec.push_back(aux_num);
                }

                if(std::isalpha(line[1]))
                {
                        switch (line[1])
                        {
                        case 'x':
                                camera.hx = std::stof(line.substr(line.find("=") + 1));
                                break;
                        case 'y':
                                camera.hy = std::stof(line.substr(line.find("=") + 1));
                                break;
                        default:
                                break;
                        }
                }

                switch(line[0])
                {
                        case 'N':
                                camera.N.x = aux_vec[0];
                                camera.N.y = aux_vec[1];
                                camera.N.z = aux_vec[2];
                                break;
                        case 'V':
                                camera.V.x = aux_vec[0];
                                camera.V.y = aux_vec[1];
                                camera.V.z = aux_vec[2];
                                break;
                        case 'C':
                                camera.C.x = aux_vec[0];
                                camera.C.y = aux_vec[1];
                                camera.C.z = aux_vec[2];
                                break;
                        case 'd':
                                camera.d = aux_vec[0];
                                break;
                        default:
                                break;
                }

                aux_vec.clear();
        }
        
}

void read_lighting_file(Lighting& lighting)
{
        std::ifstream file("lighting.txt");
        std::string line;
        while (std::getline(file, line)) 
        {
                std::istringstream iss(line.substr(line.find("=") + 1));
                std::vector <float> aux_vec;
                float aux_num;
            
                while(iss >> aux_num)
                {
                        aux_vec.push_back(aux_num);
                }

                switch(line[0])
                {
                        case 'I':
                            if(line[1] == 'l')
                            {
                                lighting.Il.x = aux_vec[0];
                                lighting.Il.y = aux_vec[1];
                                lighting.Il.z = aux_vec[2];
                            }
                            else
                            {
                                lighting.Iamb.x = aux_vec[0];
                                lighting.Iamb.y = aux_vec[1];
                                lighting.Iamb.z = aux_vec[2];
                            }
                            break;

                        case 'K':
                            switch(line[1])
                            {
                                case 'a':
                                    lighting.Ka = aux_vec[0];
                                    break;
                                case 'd':
                                    lighting.Kd.x = aux_vec[0];
                                    lighting.Kd.y = aux_vec[1];
                                    lighting.Kd.z = aux_vec[2];
                                    break;
                                case 's':
                                    lighting.Ks = aux_vec[0];
                                    break;
                                default:
                                    break;
                            }
                            break;
                        case 'O':
                            lighting.Od.x = aux_vec[0];
                            lighting.Od.y = aux_vec[1];
                            lighting.Od.z = aux_vec[2];
                            break;
                        case 'P':
                            lighting.Pl.x = aux_vec[0];
                            lighting.Pl.y = aux_vec[1];
                            lighting.Pl.z = aux_vec[2];
                        default:
                            lighting.eta = aux_vec[0];
                            break;
                }

                aux_vec.clear();
        }
        
}

void load_object(Screen& screen, Mesh& mesh, std::vector<float> numbers)
{
        //save vertices and triangle numbers
        int num_vertex = numbers[0], num_triang = numbers[1];

        //save points in points_vec
        std::vector<Vec3> points_vec;
        for (int i = 2; i < (num_vertex * 3 + 2); i += 3)
        {
                points_vec.push_back({numbers[i], numbers[i + 1], numbers[i + 2]});
        }
        
        //save triangles order in triangles_vec
        std::vector<float> triangles_vec = std::vector<float>(numbers.begin() + (num_vertex * 3 + 2), numbers.end());

        // save triangles in mesh struct according to triangles_vec
        for (int i = 0; i < triangles_vec.size(); i += 3)
        {
                mesh.triangles.push_back(
                    {{{points_vec[triangles_vec[i] - 1],
                       points_vec[triangles_vec[i + 1] - 1],
                       points_vec[triangles_vec[i + 2] - 1]}},
                       {{points_vec[triangles_vec[i] - 1],
                       points_vec[triangles_vec[i + 1] - 1],
                       points_vec[triangles_vec[i + 2] - 1]}}});
                
        }
}

void load_camera(Camera& camera)
{
        Vec3 V = camera.V, N = camera.N, C = camera.C;
        float d = camera.d, hx = camera.hx, hy = camera.hy;
        
        Vec3 v_ort = ortoghonalize(V, N);
        Vec3 U = vec_crossProduct(N, v_ort);

        camera.U = vec_normalise(U);
        camera.V = vec_normalise(v_ort);
        camera.N = vec_normalise(N);
}

//change triangle points from world coordinates to view coordinates
void change_basis(Camera camera, Mesh& mesh, std::unordered_map<std::string, std::vector<Vec3>>& hashmap)
{
        mat3x3 matrix = change_basis_matrix(camera.U, camera.V, camera.N);
        Vec3 c = camera.C;
        Vec3 p_sub0, p_sub1, p_sub2, m0, m1, m2;
        Vec3 n_sub1, n_sub2, normal;
        std::string key0, key1, key2;
        for(auto& triangle: mesh.triangles)
        {
                p_sub0 = vec_sub(triangle.points_screen[0], c);
                m0 = matrix_multVec(matrix, p_sub0);
                triangle.points_screen[0].x = m0.x;
                triangle.points_screen[0].y = m0.y;
                triangle.points_screen[0].z = m0.z;

                p_sub1 = vec_sub(triangle.points_screen[1], c);
                m1 = matrix_multVec(matrix, p_sub1);
                triangle.points_screen[1].x = m1.x;
                triangle.points_screen[1].y = m1.y;
                triangle.points_screen[1].z = m1.z;

                p_sub2 = vec_sub(triangle.points_screen[2], c);
                m2 = matrix_multVec(matrix, p_sub2);
                triangle.points_screen[2].x = m2.x;
                triangle.points_screen[2].y = m2.y;
                triangle.points_screen[2].z = m2.z;

                n_sub1 = vec_sub(triangle.points[1], triangle.points[0]);
                n_sub2 = vec_sub(triangle.points[2], triangle.points[0]);
                normal = vec_normalise(vec_crossProduct(n_sub1, n_sub2));
                triangle.normal.x = normal.x;
                triangle.normal.y = normal.y;
                triangle.normal.z = normal.z;

                key0 = std::to_string((int)triangle.points[0].x) 
                        + ',' + std::to_string((int)triangle.points[0].y) 
                        + ',' + std::to_string((int)triangle.points[0].z);
                
                key1 = std::to_string((int)triangle.points[1].x) 
                        + ',' + std::to_string((int)triangle.points[1].y) 
                        + ',' + std::to_string((int)triangle.points[1].z);
                
                key2 = std::to_string((int)triangle.points[2].x) 
                        + ',' + std::to_string((int)triangle.points[2].y) 
                        + ',' + std::to_string((int)triangle.points[2].z);

                hashmap[key0].push_back(triangle.normal);
                hashmap[key1].push_back(triangle.normal);
                hashmap[key2].push_back(triangle.normal);
        }
}

int xnorm_toScreenCoords(float x)
{
        float div = (x + 1)/2;

        return floor(div*WIDTH + 0.5);
}

int ynorm_toScreenCoords(float y)
{
        float div = (y + 1)/2;
        
        return floor(HEIGHT - div*HEIGHT + 0.5);
}

//projection function to change norm coordinates to screen coordinates
void projection(Camera camera, Mesh& mesh)
{
        float xs0, ys0, xs1, ys1, xs2, ys2;
        float d = camera.d, hx = camera.hx, hy = camera.hy;
        float d_hx = d/hx, d_hy = d/hy;
        
        for(auto& triangle: mesh.triangles)
        {
                xs0 = d_hx * (triangle.points_screen[0].x/triangle.points_screen[0].z);
                ys0 = d_hy * (triangle.points_screen[0].y/triangle.points_screen[0].z);
                
                xs1 = d_hx * (triangle.points_screen[1].x/triangle.points_screen[1].z);
                ys1 = d_hy * (triangle.points_screen[1].y/triangle.points_screen[1].z);

                xs2 = d_hx * (triangle.points_screen[2].x/triangle.points_screen[2].z);
                ys2 = d_hy * (triangle.points_screen[2].y/triangle.points_screen[2].z);

                triangle.points_screen[0].x = xnorm_toScreenCoords(xs0);
                triangle.points_screen[0].y = ynorm_toScreenCoords(ys0);

                triangle.points_screen[1].x = xnorm_toScreenCoords(xs1);
                triangle.points_screen[1].y = ynorm_toScreenCoords(ys1);

                triangle.points_screen[2].x = xnorm_toScreenCoords(xs2);
                triangle.points_screen[2].y = ynorm_toScreenCoords(ys2);
        }
}

void start_zbuffer(float (&zbuffer)[WIDTH][HEIGHT])
{
        for (int i = 0; i < WIDTH; i++)
        {
                for (int j = 0; j < HEIGHT; j++)
                {
                        zbuffer[i][j] = std::numeric_limits<float>::max();
                }
        }
}

void top_bottom(Screen& screen, Lighting& light, std::vector<Vec3> normals, float x0, float y0, float z0, float x1, float y1, float z1, float x2, float y2, float z2) {
        Vec3 Ia, Id, Is, color;
        Vec3 v, l, r, n;
        float dotProd_NL, dotProd_VR;

        float invA1 = (x1 - x0) / (y1 - y0);
        float invA2 = (x2 - x0) / (y2 - y0);

        float xMin = x0;
        float xMax = x0;

        int y = std::round(y0-1);

        xMin -= invA1;
        xMax -= invA2;

        while (y >= std::max(y1, y2)) {

                for (int i = std::round(xMin); i <= xMax; i++)
                {
                        Coeff coefficients = barycentric_coeff(i, y, x0, y0, x1, y1, x2, y2);
                        Vec3 p = {x0 * coefficients.alpha + x1 * coefficients.beta + x2 * coefficients.gamma,
                                  y0 * coefficients.alpha + y1 * coefficients.beta + y2 * coefficients.gamma,
                                  z0 * coefficients.alpha + z1 * coefficients.beta + z2 * coefficients.gamma};

                        if (p.z > 0 && p.z < zbuffer[i][y])
                        {
                                n = vec_normalise({normals[0].x*coefficients.alpha + normals[1].x*coefficients.beta + normals[2].x*coefficients.gamma,
                                                   normals[0].y*coefficients.alpha + normals[1].y*coefficients.beta + normals[2].y*coefficients.gamma,
                                                   normals[0].z*coefficients.alpha + normals[1].z*coefficients.beta + normals[2].z*coefficients.gamma});

                                v = vec_normalise({-p.x, -p.y, -p.z});
                                l = vec_normalise(vec_sub(light.Pl, p));
                                dotProd_NL = vec_dotProduct(n, l);
                                r = vec_sub(vec_mul(n, (2 * dotProd_NL)), l);
                                dotProd_VR = vec_dotProduct(v, r);
                                
                                Ia = vec_mul(light.Iamb, light.Ka);
                                Id = vec_mul(component_mul(light.Kd, light.Od, light.Il), dotProd_NL);
                                Is = vec_mul(light.Il, (std::pow(dotProd_VR, light.eta) * light.Ks));
                                
                                if ((dotProd_NL < 0))
                                {
                                        if (vec_dotProduct(v, n) < 0)
                                        {
                                                n.x *= -1;
                                                n.y *= -1;
                                                n.z *= -1;
                                        }
                                        else
                                        {
                                                Id = Is = {0, 0, 0};
                                        }
                                }

                                if (dotProd_VR < 0)
                                {
                                        Is = {0, 0, 0};
                                }

                                color = {Ia.x + Id.x + Is.x,
                                         Ia.y + Id.y + Is.y,
                                         Ia.z + Id.z + Is.z};

                                color = {color.x > 255 ? 255 : color.x,
                                         color.y > 255 ? 255 : color.y,
                                         color.z > 255 ? 255 : color.z};

                                zbuffer[i][y] = p.z;
                                screen.pixel(i, y, color.x, color.y, color.z);
                        }
                }

                y--;
                xMin -= invA1;
                xMax -= invA2;
        }
}

void bottom_up(Screen& screen,Lighting& light, std::vector<Vec3> normals, float x0, float y0, float z0, float x1, float y1, float z1, float x2, float y2, float z2)
{
        Vec3 Ia, Id, Is, color;
        Vec3 v, l, r, n;
        float dotProd_NL, dotProd_VR;
        float area = 0.5 * (x0 * (y1 - y2) + x1 * (y2 - y0) + x2 * (y0 - y1));

        if (area == 0) {
                return;
        }

        float invA1 = (x1 - x0) / (y1 - y0);
        float invA2 = (x2 - x0) / (y2 - y0);

        if (invA1 > invA2) {
                normals = {normals[0], normals[2], normals[1]};
                bottom_up(screen, light, normals, x0, y0, z0, x2, y2, z2, x1, y1, z1);
                return;
        }

        float xMin = x0;
        float xMax = xMin;

        int y = y0+1;

        xMin += invA1;
        xMax += invA2;

        while (y <= std::min(y1, y2)) {
                
                for (int i = (int)std::round(xMin); i <= xMax; i++)
                {
                        Coeff coefficients = barycentric_coeff(i, y, x0, y0, x1, y1, x2, y2);
                        Vec3 p = {x0*coefficients.alpha + x1*coefficients.beta + x2*coefficients.gamma,
                                  y0*coefficients.alpha + y1*coefficients.beta + y2*coefficients.gamma,
                                  z0*coefficients.alpha + z1*coefficients.beta + z2*coefficients.gamma};

                        if(p.z > 0 && p.z < zbuffer[i][y])
                        {
                                n = vec_normalise({normals[0].x*coefficients.alpha + normals[1].x*coefficients.beta + normals[2].x*coefficients.gamma,
                                                   normals[0].y*coefficients.alpha + normals[1].y*coefficients.beta + normals[2].y*coefficients.gamma,
                                                   normals[0].z*coefficients.alpha + normals[1].z*coefficients.beta + normals[2].z*coefficients.gamma});

                                v = vec_normalise({-p.x, -p.y, -p.z});
                                l = vec_normalise(vec_sub(light.Pl, p));
                                dotProd_NL = vec_dotProduct(n,l);
                                r = vec_sub(vec_mul(n,(2*dotProd_NL)),l);
                                dotProd_VR = vec_dotProduct(v,r);
                                
                                Ia = vec_mul(light.Iamb, light.Ka);
                                Id = vec_mul(component_mul(light.Kd, light.Od, light.Il), dotProd_NL);
                                Is = vec_mul(light.Il,(std::pow(dotProd_VR, light.eta)*light.Ks));
                                
                                if((dotProd_NL < 0))
                                {
                                        if(vec_dotProduct(v,n) < 0)
                                        {
                                                n.x *= -1;
                                                n.y *= -1;
                                                n.z *= -1;   
                                        }
                                        else
                                        {
                                                Id = {0,0,0};
                                                Is = {0,0,0};
                                        }  
                                }

                                if(dotProd_VR < 0)
                                {
                                        Is = {0,0,0};
                                }

                                color = {Ia.x + Id.x + Is.x,
                                         Ia.y + Id.y + Is.y,
                                         Ia.z + Id.z + Is.z};

                                color = {color.x > 255 ? 255 : color.x,
                                         color.y > 255 ? 255 : color.y,
                                         color.z > 255 ? 255 : color.z};

                                zbuffer[i][y] = p.z;
                                screen.pixel(i, y, color.x, color.y, color.z);
                        }
                        
                }

                y++;
                xMin += invA1;
                xMax += invA2;
        }

        if (y <= std::max(y1, y2)) {

                if (y1 > y2) {
                        top_bottom(screen, light, normals, x1, y1, z1, (int) std::round(xMin), y, z0, x2, y2, z2);
                } else {
                        top_bottom(screen, light, normals, x2, y2, z2, x1, y1, z1, (int) std::round(xMax), y, z0);
                }
        }
}


void scanline(Screen &screen, Mesh mesh, Lighting& light, std::unordered_map<std::string, std::vector<Vec3>> hashmap)
{
        std::vector<Vec3> normals;
        Vec3 n0, n1, n2;
        std::string key0, key1, key2;
        float sum_x, sum_y, sum_z;

        for(auto& triangle: mesh.triangles)
        {

                std::sort(triangle.points_screen.begin(), triangle.points_screen.end(), [](const Vec3& a, const Vec3& b) {
                        return a.y < b.y;
                });
                
                std::sort(triangle.points.begin(), triangle.points.end(), [](const Vec3& a, const Vec3& b) {
                        return a.y < b.y;
                });
        
                key0 = std::to_string((int)triangle.points[0].x) 
                        + ',' + std::to_string((int)triangle.points[0].y) 
                        + ',' + std::to_string((int)triangle.points[0].z);
                
                key1 = std::to_string((int)triangle.points[1].x) 
                        + ',' + std::to_string((int)triangle.points[1].y) 
                        + ',' + std::to_string((int)triangle.points[1].z);
                
                key2 = std::to_string((int)triangle.points[2].x) 
                        + ',' + std::to_string((int)triangle.points[2].y) 
                        + ',' + std::to_string((int)triangle.points[2].z);
                
                for(auto& n: hashmap[key0])
                {
                        sum_x += n.x;
                        sum_y += n.y;
                        sum_z += n.z;
                }

                sum_x /= hashmap[key0].size();
                sum_y /= hashmap[key0].size();
                sum_z /= hashmap[key0].size();
                n0 = vec_normalise({sum_x, sum_y, sum_z}); 

                sum_x = sum_y = sum_z = 0;

                for(auto& n: hashmap[key1])
                {
                        sum_x += n.x;
                        sum_y += n.y;
                        sum_z += n.z;
                }

                sum_x /= hashmap[key1].size();
                sum_y /= hashmap[key1].size();
                sum_z /= hashmap[key1].size();
                n1 = vec_normalise({sum_x, sum_y, sum_z});

                sum_x = sum_y = sum_z = 0;

                for(auto& n: hashmap[key2])
                {
                        sum_x += n.x;
                        sum_y += n.y;
                        sum_z += n.z;
                }

                sum_x /= hashmap[key2].size();
                sum_y /= hashmap[key2].size();
                sum_z /= hashmap[key2].size();
                n2 = vec_normalise({sum_x, sum_y, sum_z}); 

                sum_x = sum_y = sum_z = 0;

                normals.push_back(n0);
                normals.push_back(n1);
                normals.push_back(n2); 

                bottom_up(screen, light, normals,
                          triangle.points_screen[0].x, triangle.points_screen[0].y, triangle.points_screen[0].z,
                          triangle.points_screen[1].x, triangle.points_screen[1].y, triangle.points_screen[1].z,
                          triangle.points_screen[2].x, triangle.points_screen[2].y, triangle.points_screen[2].z);

                normals.clear();
        }

}

void input(Screen& screen, Mesh& mesh, Camera& camera, Lighting& light, std::unordered_map<std::string, std::vector<Vec3>> hashmap, std::vector <float>& numbers)
{
        SDL_Event ev;
        
        while(SDL_PollEvent(&ev))
        {
                if(ev.type == SDL_KEYDOWN && ev.key.keysym.sym == SDLK_r)
                {
                        numbers.clear();
                        screen.clear();
                        mesh = {};

                        read_file(numbers); 
                        load_object(screen, mesh, numbers);
                        start_zbuffer(zbuffer);
                        read_camera_file(camera);
                        load_camera(camera);
                        change_basis(camera, mesh, hashmap);
                        projection(camera, mesh);
                        scanline(screen, mesh, light, hashmap);
                }
                else if (ev.type == SDL_QUIT)
                {
                        screen.quit();
                }
                
        }

}

int main(int argc, char *argv[])
{       
        Screen screen;
        Mesh mesh;
        Camera camera;
        Lighting light;
        std::vector<float> numbers;     
        std::unordered_map<std::string, std::vector<Vec3>> hashmap; 

        read_file(numbers);
        load_object(screen, mesh, numbers);
        start_zbuffer(zbuffer);
        read_camera_file(camera);
        read_lighting_file(light);
        load_camera(camera);
        change_basis(camera, mesh, hashmap);
        projection(camera, mesh);
        scanline(screen, mesh, light, hashmap);
        
        while (true)
        {
                screen.show();
                input(screen, mesh, camera, light, hashmap, numbers);
        }
        
        return 0;
}