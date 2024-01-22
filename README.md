# Projeto de Computação Gráfica

Este é o projeto da cadeira de computação gráfica básica. A linguagem escolhida foi C++ e a biblioteca externa utilizada foi SDL2, para usar a função de pintar um pixel na tela. O projeto consiste em renderizar um objeto na tela a partir de um arquivo que descreve uma malha com todos os pontos e triângulos formados pelos pontos. Para a renderização dos triângulos é utilizado o algoritmo de [Scan Line](https://en.wikipedia.org/wiki/Scanline_rendering) junto ao [Z-Buffer](https://en.wikipedia.org/wiki/Z-buffering), e para fazer a iluminação e tonalização é utilizado o algoritmo de [Phong](https://www.scratchapixel.com/lessons/3d-basic-rendering/phong-shader-BRDF/phong-illumination-models-brdf.html).

## Compilação
O código foi compilado em um sistema Windows com o comando abaixo:
```console
$ g++ -I src/include -I src/headers -L src/lib -o projectfinal.exe main.cpp -lmingw32 -lSDL2main -lSDL2
```
Para compilar em outro sistema (Linux, MacOS, etc) podem ser necessárias flags diferentes.

## Execução
O executável **projectfinal.exe** está lendo as informações dos arquivos **object.byu** (o objeto guardado nele é um vaso) e **camera.txt**. Para atualizar o objeto ou a câmera, após alterar o respectivo arquivo, aperte a tecla **r** do teclado uma vez.