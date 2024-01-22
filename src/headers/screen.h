#include <SDL2/SDL.h>
#include <vector>
#include <algorithm>

#define WIDTH 1280
#define HEIGHT 960

class Screen
{
        SDL_Event e;
        SDL_Event reload;
        SDL_Window* window;
        SDL_Renderer* renderer;

        public:
        Screen()
        {
                SDL_Init(SDL_INIT_VIDEO);
                SDL_CreateWindowAndRenderer(
                    WIDTH, HEIGHT, 0, &window, &renderer);
                SDL_RenderSetScale(renderer, 1, 1);
                SDL_SetRenderDrawColor(renderer, 0,0,0,255);
                SDL_RenderClear(renderer);
        }

        void pixel(int x, int y, int r, int g, int b)
        {
                SDL_SetRenderDrawColor(renderer, r, g, b, 255);
                SDL_RenderDrawPointF(renderer, x, y); 
        }

        void show()
        {
                SDL_RenderPresent(renderer);
        }

        void clear()
        {
                SDL_SetRenderDrawColor(renderer,0,0,0,255);
                SDL_RenderClear(renderer);
                SDL_RenderPresent(renderer);
        }

        void quit()
        {
                SDL_Quit();
                exit(0);
        }
};