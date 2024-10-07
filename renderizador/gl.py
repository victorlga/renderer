#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

# pylint: disable=invalid-name

"""
Biblioteca Gráfica / Graphics Library.

Desenvolvido por: Victor Luís Gama de Assis
Disciplina: Computação Gráfica
Data: Aug 13th, 2024
"""

import time         # Para operações com tempo
import gpu          # Simula os recursos de uma GPU
import math         # Funções matemáticas
import numpy as np  # Biblioteca do Numpy

class GL:
    """Classe que representa a biblioteca gráfica (Graphics Library)."""

    width = 800   # largura da tela
    height = 600  # altura da tela
    near = 0.01   # plano de corte próximo
    far = 1000    # plano de corte distante

    @staticmethod
    def setup(width, height, near=0.01, far=1000):
        """Definr parametros para câmera de razão de aspecto, plano próximo e distante."""
        GL.width = width * 2
        GL.height = height * 2
        GL.near = near
        GL.far = far
        GL.perspective_matrix = None
        GL.transformation_stack = []

    @staticmethod
    def polypoint2D(point, colors):
        """Função usada para renderizar Polypoint2D."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry2D.html#Polypoint2D
        # Nessa função você receberá pontos no parâmetro point, esses pontos são uma lista
        # de pontos x, y sempre na ordem. Assim point[0] é o valor da coordenada x do
        # primeiro ponto, point[1] o valor y do primeiro ponto. Já point[2] é a
        # coordenada x do segundo ponto e assim por diante. Assuma a quantidade de pontos
        # pelo tamanho da lista e assuma que sempre vira uma quantidade par de valores.
        # O parâmetro colors é um dicionário com os tipos cores possíveis, para o Polypoint2D
        # você pode assumir inicialmente o desenho dos pontos com a cor emissiva (emissiveColor).

        i = 0
        color = list(map(lambda x: round(x * 255), colors['emissiveColor']))  # Convertendo a cor para valores entre 0 e 255
        while i < len(point):
            gpu.GPU.draw_pixel([int(point[i]), int(point[i+1])], gpu.GPU.RGB8, color)  # Desenhando cada ponto
            i += 2

    @staticmethod
    def polyline2D(lineSegments, colors):
        """Função usada para renderizar Polyline2D."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry2D.html#Polyline2D
        # Nessa função você receberá os pontos de uma linha no parâmetro lineSegments, esses
        # pontos são uma lista de pontos x, y sempre na ordem. Assim point[0] é o valor da
        # coordenada x do primeiro ponto, point[1] o valor y do primeiro ponto. Já point[2] é
        # a coordenada x do segundo ponto e assim por diante. Assuma a quantidade de pontos
        # pelo tamanho da lista. A quantidade mínima de pontos são 2 (4 valores), porém a
        # função pode receber mais pontos para desenhar vários segmentos. Assuma que sempre
        # vira uma quantidade par de valores.
        # O parâmetro colors é um dicionário com os tipos cores possíveis, para o Polyline2D
        # você pode assumir inicialmente o desenho das linhas com a cor emissiva (emissiveColor).

        print("Polyline2D : lineSegments = {0}".format(lineSegments)) # imprime no terminal
        print("Polyline2D : colors = {0}".format(colors)) # imprime no terminal as cores
        
        x0, y0, x1, y1 = map(int, lineSegments)  # Convertendo as coordenadas para inteiros
        color = list(map(lambda x: round(x * 255), colors['emissiveColor']))  # Convertendo a cor para valores entre 0 e 255

        if y0 == y1:
            # Caso especial de linha horizontal
            for u in range(x0, x1+1):
                gpu.GPU.draw_pixel([u, y0], gpu.GPU.RGB8, color)  # Desenha a linha horizontal
            return
    
        if y0 > y1:
            # Se y0 > y1, inverta os pontos para garantir que y0 esteja sempre abaixo de y1
            y0, y1 = y1, y0
            x0, x1 = x1, x0

        s = (x1 - x0) / (y1 - y0)  # Calcula a inclinação da linha

        if abs(s) < 1:
            # Se a inclinação for menor que 1, iteramos sobre y e calculamos x
            u = x0
            for v in range(y0, y1+1):
                gpu.GPU.draw_pixel([round(u), v], gpu.GPU.RGB8, color)  # Desenha o ponto calculado
                u += s
            return

        # Caso contrário, iteramos sobre x e calculamos y
        y0, y1 = y1, y0
        x0, x1 = x1, x0
    
        s **= -1  # Inverte a inclinação para a iteração em x
        v = y0
        for u in range(x0, x1+1):
            gpu.GPU.draw_pixel([u, round(v)], gpu.GPU.RGB8, color)  # Desenha o ponto calculado
            v += s

    @staticmethod
    def circle2D(radius, colors):
        """Função usada para renderizar Circle2D."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry2D.html#Circle2D
        # Nessa função você receberá um valor de raio e deverá desenhar o contorno de
        # um círculo.
        # O parâmetro colors é um dicionário com os tipos cores possíveis, para o Circle2D
        # você pode assumir o desenho das linhas com a cor emissiva (emissiveColor).

        # Convertendo a cor para valores entre 0 e 255
        color = (np.array(colors['emissiveColor']) * 255).astype(int)

        # Cria um array de ângulos de 0 a 2*pi
        angles = np.linspace(0, 2 * np.pi, num=360, endpoint=False)

        # Calcula as coordenadas x e y dos pontos na circunferência
        x_points = np.round(np.cos(angles) * radius).astype(int)
        y_points = np.round(np.sin(angles) * radius).astype(int)

        # Combina as coordenadas em pares de pontos
        points = np.vstack((x_points, y_points)).T

        # Filtra os pontos para garantir que estão dentro do framebuffer
        valid_points = points[(points[:, 0] >= 0) & (points[:, 0] < GL.width) & 
                            (points[:, 1] >= 0) & (points[:, 1] < GL.height)]

        # Desenha cada ponto na tela
        for point in valid_points:
            gpu.GPU.draw_pixel(point.tolist(), gpu.GPU.RGB8, color.tolist())
    
    
    @staticmethod
    def is_inside(vertices, point):
        # Utiliza numpy para cálculos vetorizados do produto vetorial
        x0, y0, x1, y1, x2, y2 = np.array(vertices).flatten()

        v0 = np.array([x2 - x0, y2 - y0])
        v1 = np.array([x1 - x0, y1 - y0])
        v2 = np.array(point) - np.array([x0, y0])

        dot00 = np.dot(v0, v0)
        dot01 = np.dot(v0, v1)
        dot02 = np.dot(v0, v2)
        dot11 = np.dot(v1, v1)
        dot12 = np.dot(v1, v2)

        invDenom = 1 / (dot00 * dot11 - dot01 * dot01)
        u = (dot11 * dot02 - dot01 * dot12) * invDenom
        v = (dot00 * dot12 - dot01 * dot02) * invDenom

        return (u >= 0) & (v >= 0) & (u + v <= 1)

    
    @staticmethod
    def compute_barycentric_coordinates(vertices, point):
        """Calcula os pesos baricêntricos (α, β, γ) para um ponto em relação a um triângulo usando as fórmulas da imagem."""

        # Extraindo as coordenadas dos vértices
        xA, yA = vertices[0]
        xB, yB = vertices[1]
        xC, yC = vertices[2]

        # Extraindo a coordenada do ponto
        x, y = point

        # Cálculo do denominador (determinante)
        alpha_denom = -(xA - xB) * (yC - yB) + (yA - yB) * (xC - xB)
        beta_denom = -(xB - xC) * (yA - yC) + (yB - yC) * (xA - xC)

        # Cálculo de alpha
        alpha = (-(x - xB) * (yC - yB) + (y - yB) * (xC - xB)) / alpha_denom

        # Cálculo de beta
        beta = (-(x - xC) * (yA - yC) + (y - yC) * (xA - xC)) / beta_denom

        # Cálculo de gamma
        gamma = 1.0 - alpha - beta

        return alpha, beta, gamma


    @staticmethod
    def rasterize_triangle(vertices, colors, z_coords=None, transparency=None):
        
        """Função que rasteriza qualquer triângulo."""
        # Calcula o bounding box do triângulo
        min_x, min_y = np.min(vertices, axis=0).astype(int)
        max_x, max_y = np.max(vertices, axis=0).astype(int)

        # Itera dentro do bounding box e desenha os pixels
        for x in range(min_x, max_x + 1):
            for y in range(min_y, max_y + 1):
                if 0 <= x < GL.width and 0 <= y < GL.height:
                    if not GL.is_inside(vertices, [x, y]):
                        continue
                    
                    if z_coords:
                        alpha, beta, gamma = GL.compute_barycentric_coordinates(vertices, (x, y))

                        interpolated_color = (
                            alpha * (colors[0] / z_coords[0][0]) +
                            beta  * (colors[1] / z_coords[1][0]) +
                            gamma * (colors[2] / z_coords[2][0])
                        )


                        z_interpolated = 1 / (alpha / z_coords[0][0] + beta / z_coords[1][0] + gamma / z_coords[2][0])
                        interpolated_color *= z_interpolated

                        # Converte a cor interpolada para inteiro
                        interpolated_color = interpolated_color.astype(int)
                        

                        new_depht = 1 / (alpha / z_coords[0][1] + beta / z_coords[1][1] + gamma / z_coords[2][1])
                        if gpu.GPU.read_pixel([x, y], gpu.GPU.DEPTH_COMPONENT32F) > new_depht:
                            gpu.GPU.draw_pixel([x, y], gpu.GPU.DEPTH_COMPONENT32F, [new_depht])
                            
                            if transparency:
                                current_color = gpu.GPU.read_pixel([x, y], gpu.GPU.RGB8) * transparency
                                new_color = interpolated_color * (1 - transparency)
                                interpolated_color = current_color + new_color
                                
                            gpu.GPU.draw_pixel([x, y], gpu.GPU.RGB8, interpolated_color.tolist())

                    else:
                        gpu.GPU.draw_pixel([x, y], gpu.GPU.RGB8, colors.tolist())

    @staticmethod
    def triangleSet2D(vertices, colors):
        """Função usada para renderizar TriangleSet2D."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry2D.html#TriangleSet2D
        # Nessa função você receberá os vertices de um triângulo no parâmetro vertices,
        # esses pontos são uma lista de pontos x, y sempre na ordem. Assim point[0] é o
        # valor da coordenada x do primeiro ponto, point[1] o valor y do primeiro ponto.
        # Já point[2] é a coordenada x do segundo ponto e assim por diante. Assuma que a
        # quantidade de pontos é sempre multiplo de 3, ou seja, 6 valores ou 12 valores, etc.
        # O parâmetro colors é um dicionário com os tipos cores possíveis, para o TriangleSet2D
        # você pode assumir inicialmente o desenho das linhas com a cor emissiva (emissiveColor).

        colors = np.array(colors['emissiveColor']) * 255  # Convertendo a cor para valores entre 0 e 255
        colors = colors.astype(int)
        vertices = np.array(vertices).reshape(3, 2) * 2

        GL.rasterize_triangle(vertices, colors)
    
    @staticmethod
    def project_vertex(vertex):
        """Aplica a transformação e projeção a um vértice."""
        vertex_homogeneous = np.append(vertex, 1)  # Coordenadas homogêneas

        model_matrix = GL.transformation_stack[-1] if GL.transformation_stack else np.identity(4)
        transformed_vertex = model_matrix @ vertex_homogeneous
        projected_vertex = GL.perspective_matrix @ transformed_vertex
        z_camera_space = projected_vertex[2]

        # Realizamos a divisão por w
        projected_vertex /= projected_vertex[3]

        # Convertendo para coordenadas de tela
        x = (projected_vertex[0] + 1) * 0.5 * GL.width
        y = (1 - (projected_vertex[1] + 1) * 0.5) * GL.height
        return [x, y], [z_camera_space, projected_vertex[2]]

    @staticmethod
    def triangleSet(point, colors):
        """Função usada para renderizar TriangleSet."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/rendering.html#TriangleSet
        # Nessa função você receberá pontos no parâmetro point, esses pontos são uma lista
        # de pontos x, y, e z sempre na ordem. Assim point[0] é o valor da coordenada x do
        # primeiro ponto, point[1] o valor y do primeiro ponto, point[2] o valor z da
        # coordenada z do primeiro ponto. Já point[3] é a coordenada x do segundo ponto e
        # assim por diante.
        # No TriangleSet os triângulos são informados individualmente, assim os três
        # primeiros pontos definem um triângulo, os três próximos pontos definem um novo
        # triângulo, e assim por diante.
        # O parâmetro colors é um dicionário com os tipos cores possíveis, você pode assumir
        # inicialmente, para o TriangleSet, o desenho das linhas com a cor emissiva
        # (emissiveColor), conforme implementar novos materias você deverá suportar outros
        # tipos de cores.
        
        emissiveColor = np.array(colors['emissiveColor']) * 255  # Convertendo a cor para valores entre 0 e 255
        emissiveColor = emissiveColor.astype(int)
        emissiveColor = [emissiveColor] * 3
        transparency = colors['transparency']

        num_points = len(point) // 3
        points = np.array(point).reshape(num_points, 3)
        
        for i in range(0, num_points, 3):
            p1, z1 = GL.project_vertex(points[i])
            p2, z2 = GL.project_vertex(points[i+1])
            p3, z3 = GL.project_vertex(points[i+2])
            vertices = np.array(p1+p2+p3).reshape(3, 2)

            # Utiliza a função rasterize_triangle para desenhar o triângulo
            GL.rasterize_triangle(vertices, emissiveColor, z_coords=[z1, z2, z3], transparency=transparency)
    
    @staticmethod
    def quaternion_to_matrix(q):
        """Converte um quaternion em uma matriz de rotação 4x4."""
        w, x, y, z = q
        return np.array([
            [1 - 2*y*y - 2*z*z, 2*x*y - 2*z*w, 2*x*z + 2*y*w, 0],
            [2*x*y + 2*z*w, 1 - 2*x*x - 2*z*z, 2*y*z - 2*x*w, 0],
            [2*x*z - 2*y*w, 2*y*z + 2*x*w, 1 - 2*x*x - 2*y*y, 0],
            [0, 0, 0, 1]
        ])

    @staticmethod
    def viewpoint(position, orientation, fieldOfView):
        """Função usada para renderizar (na verdade coletar os dados) de Viewpoint."""
        # Na função de viewpoint você receberá a posição, orientação e campo de visão da
        # câmera virtual. Use esses dados para poder calcular e criar a matriz de projeção
        # perspectiva para poder aplicar nos pontos dos objetos geométricos.

        aspect_ratio = GL.width / GL.height
        near = GL.near
        far = GL.far

        # Matriz de projeção perspectiva
        top = near * np.tan(fieldOfView / 2)
        right = top * aspect_ratio

        perspective_matrix = np.array([
            [near/right, 0, 0, 0],
            [0, near/top, 0, 0],
            [0, 0, -(far+near)/(far-near), -2*far*near/(far-near)],
            [0, 0, -1, 0]
        ])

        # Matriz de rotação da câmera a partir de `orientation`
        half_angle = orientation[3] / 2
        x, y, z = orientation[:3]
        c, s = np.cos(half_angle), np.sin(half_angle)
        q = np.array([
            c,
            s * x,
            s * y,
            s * z
        ])
        rotation_matrix = GL.quaternion_to_matrix(q)

        # Matriz de translação da câmera a partir de `position`
        translation_matrix = np.linalg.inv(np.array([
            [1, 0, 0, position[0]],
            [0, 1, 0, position[1]],
            [0, 0, 1, position[2]],
            [0, 0, 0, 1]
        ]))

        # Combinando as transformações
        GL.perspective_matrix = perspective_matrix @ rotation_matrix @ translation_matrix


    @staticmethod
    def transform_in(translation, scale, rotation):
        """Função usada para renderizar (na verdade coletar os dados) de Transform."""
        # A função transform_in será chamada quando se entrar em um nó X3D do tipo Transform
        # do grafo de cena. Os valores passados são a escala em um vetor [x, y, z]
        # indicando a escala em cada direção, a translação [x, y, z] nas respectivas
        # coordenadas e finalmente a rotação por [x, y, z, t] sendo definida pela rotação
        # do objeto ao redor do eixo x, y, z por t radianos, seguindo a regra da mão direita.
        # Quando se entrar em um nó transform se deverá salvar a matriz de transformação dos
        # modelos do mundo para depois potencialmente usar em outras chamadas. 
        # Quando começar a usar Transforms dentre de outros Transforms, mais a frente no curso
        # Você precisará usar alguma estrutura de dados pilha para organizar as matrizes.

        current_matrix = GL.transformation_stack[-1] if GL.transformation_stack else np.identity(4)

        # Matriz de translação
        translation_matrix = np.array([
            [1, 0, 0, translation[0]],
            [0, 1, 0, translation[1]],
            [0, 0, 1, translation[2]],
            [0, 0, 0, 1]
        ])

        # Matriz de escala
        scale_matrix = np.diag([scale[0], scale[1], scale[2], 1])

        # Conversão de rotação para quaternion
        angle = rotation[3]
        axis = np.array(rotation[:3])
        axis = axis / np.linalg.norm(axis)  # Normaliza o eixo de rotação
        half_angle = angle / 2.0
        sin_half_angle = np.sin(half_angle)

        q = np.array([
            np.cos(half_angle),
            sin_half_angle * axis[0],
            sin_half_angle * axis[1],
            sin_half_angle * axis[2]
        ])

        # Converte o quaternion em uma matriz de rotação
        rotation_matrix = GL.quaternion_to_matrix(q)

        transformation_matrix = current_matrix @ translation_matrix @ rotation_matrix @ scale_matrix

        GL.transformation_stack.append(transformation_matrix)

    @staticmethod
    def transform_out():
        """Função usada para renderizar (na verdade coletar os dados) de Transform."""
        # A função transform_out será chamada quando se sair em um nó X3D do tipo Transform do
        # grafo de cena. Não são passados valores, porém quando se sai de um nó transform se
        # deverá recuperar a matriz de transformação dos modelos do mundo da estrutura de
        # pilha implementada.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        if GL.transformation_stack:
            GL.transformation_stack.pop()

    @staticmethod
    def triangleStripSet(point, stripCount, colors):
        """Função usada para renderizar TriangleStripSet."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/rendering.html#TriangleStripSet
        # A função triangleStripSet é usada para desenhar tiras de triângulos interconectados,
        # você receberá as coordenadas dos pontos no parâmetro point, esses pontos são uma
        # lista de pontos x, y, e z sempre na ordem. Assim point[0] é o valor da coordenada x
        # do primeiro ponto, point[1] o valor y do primeiro ponto, point[2] o valor z da
        # coordenada z do primeiro ponto. Já point[3] é a coordenada x do segundo ponto e assim
        # por diante. No TriangleStripSet a quantidade de vértices a serem usados é informado
        # em uma lista chamada stripCount (perceba que é uma lista). Ligue os vértices na ordem,
        # primeiro triângulo será com os vértices 0, 1 e 2, depois serão os vértices 1, 2 e 3,
        # depois 2, 3 e 4, e assim por diante. Cuidado com a orientação dos vértices, ou seja,
        # todos no sentido horário ou todos no sentido anti-horário, conforme especificado.

        colors = np.array(colors['emissiveColor']) * 255
        
        colors = colors.astype(int)
        colors = [colors] * 3

        offset = 0

        for count in stripCount:
            # Vamos percorrer os vértices da tira de triângulos
            for i in range(count - 2):
                # Pegamos três vértices consecutivos
                p1 = point[(offset + i) * 3:(offset + i) * 3 + 3]
                p2 = point[(offset + i + 1) * 3:(offset + i + 1) * 3 + 3]
                p3 = point[(offset + i + 2) * 3:(offset + i + 2) * 3 + 3]

                # Se o índice é ímpar, invertemos a ordem dos vértices para garantir a orientação correta
                if i % 2 == 1:
                    p1, p3 = p3, p1
                
                # Projeta os vértices para o espaço 2D
                p1_2d, z1 = GL.project_vertex(np.array(p1))
                p2_2d, z2 = GL.project_vertex(np.array(p2))
                p3_2d, z3 = GL.project_vertex(np.array(p3))

                # Desenhe o triângulo
                GL.rasterize_triangle([p1_2d, p2_2d, p3_2d], colors, z_coords=[z1, z2, z3])

            # Atualiza o deslocamento para a próxima tira de triângulos
            offset += count

    @staticmethod
    def indexedTriangleStripSet(point, index, colors):
        """Função usada para renderizar IndexedTriangleStripSet."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/rendering.html#IndexedTriangleStripSet
        # A função indexedTriangleStripSet é usada para desenhar tiras de triângulos
        # interconectados, você receberá as coordenadas dos pontos no parâmetro point, esses
        # pontos são uma lista de pontos x, y, e z sempre na ordem. Assim point[0] é o valor
        # da coordenada x do primeiro ponto, point[1] o valor y do primeiro ponto, point[2]
        # o valor z da coordenada z do primeiro ponto. Já point[3] é a coordenada x do
        # segundo ponto e assim por diante. No IndexedTriangleStripSet uma lista informando
        # como conectar os vértices é informada em index, o valor -1 indica que a lista
        # acabou. A ordem de conexão será de 3 em 3 pulando um índice. Por exemplo: o
        # primeiro triângulo será com os vértices 0, 1 e 2, depois serão os vértices 1, 2 e 3,
        # depois 2, 3 e 4, e assim por diante. Cuidado com a orientação dos vértices, ou seja,
        # todos no sentido horário ou todos no sentido anti-horário, conforme especificado.
        # Converte as cores para valores entre 0 e 255
        colors = np.array(colors['emissiveColor']) * 255
        colors = colors.astype(int)
        colors = [colors] * 3

        # Precisamos percorrer a lista de índices e montar os triângulos
        i = 0  # Índice para percorrer o array de índices
        
        while i < len(index) - 2:
            # Verifica se encontramos o delimitador -1
            if index[i] == -1 or index[i + 1] == -1 or index[i + 2] == -1:
                i += 1
                continue
            
            # Obtemos os três índices para o triângulo
            idx1 = index[i]
            idx2 = index[i + 1]
            idx3 = index[i + 2]

            # Obtemos os vértices correspondentes usando os índices
            p1 = point[idx1 * 3:idx1 * 3 + 3]
            p2 = point[idx2 * 3:idx2 * 3 + 3]
            p3 = point[idx3 * 3:idx3 * 3 + 3]

            # Para alternar a orientação, invertemos os vértices de triângulos ímpares
            if i % 2 == 1:
                p1, p3 = p3, p1

            # Projeta os vértices no espaço 2D
            p1_2d, z1 = GL.project_vertex(np.array(p1))
            p2_2d, z2 = GL.project_vertex(np.array(p2))
            p3_2d, z3 = GL.project_vertex(np.array(p3))
            
            # Desenha o triângulo
            GL.rasterize_triangle([p1_2d, p2_2d, p3_2d], colors, z_coords=[z1, z2, z3])

            # Incrementa o índice para a próxima sequência
            i += 1

    @staticmethod
    def indexedFaceSet(coord, coordIndex, colorPerVertex, color, colorIndex,
                       texCoord, texCoordIndex, colors, current_texture):
        """Função usada para renderizar IndexedFaceSet."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry3D.html#IndexedFaceSet
        # A função indexedFaceSet é usada para desenhar malhas de triângulos. Ela funciona de
        # forma muito simular a IndexedTriangleStripSet porém com mais recursos.
        # Você receberá as coordenadas dos pontos no parâmetro cord, esses
        # pontos são uma lista de pontos x, y, e z sempre na ordem. Assim coord[0] é o valor
        # da coordenada x do primeiro ponto, coord[1] o valor y do primeiro ponto, coord[2]
        # o valor z da coordenada z do primeiro ponto. Já coord[3] é a coordenada x do
        # segundo ponto e assim por diante. No IndexedFaceSet uma lista de vértices é informada
        # em coordIndex, o valor -1 indica que a lista acabou.
        # A ordem de conexão não possui uma ordem oficial, mas em geral se o primeiro ponto com os dois
        # seguintes e depois este mesmo primeiro ponto com o terçeiro e quarto ponto. Por exemplo: numa
        # sequencia 0, 1, 2, 3, 4, -1 o primeiro triângulo será com os vértices 0, 1 e 2, depois serão
        # os vértices 0, 2 e 3, e depois 0, 3 e 4, e assim por diante, até chegar no final da lista.
        # Adicionalmente essa implementação do IndexedFace aceita cores por vértices, assim
        # se a flag colorPerVertex estiver habilitada, os vértices também possuirão cores
        # que servem para definir a cor interna dos poligonos, para isso faça um cálculo
        # baricêntrico de que cor deverá ter aquela posição. Da mesma forma se pode definir uma
        # textura para o poligono, para isso, use as coordenadas de textura e depois aplique a
        # cor da textura conforme a posição do mapeamento. Dentro da classe GPU já está
        # implementadado um método para a leitura de imagens.

        # Conversão da cor emissiva
        emissiveColor = np.array(colors['emissiveColor']) * 255
        emissiveColor = emissiveColor.astype(int)

        i = 0
        face_indices = []
        color_indices = []

        while i < len(coordIndex):
            if coordIndex[i] == -1:
                # Face completa, processar
                if len(face_indices) >= 3:
                    # Triangularizar a face
                    for j in range(1, len(face_indices) - 1):
                        idx1 = face_indices[0]
                        idx2 = face_indices[j]
                        idx3 = face_indices[j + 1]

                        # Obter os vértices correspondentes
                        p1 = np.array(coord[idx1 * 3:idx1 * 3 + 3])
                        p2 = np.array(coord[idx2 * 3:idx2 * 3 + 3])
                        p3 = np.array(coord[idx3 * 3:idx3 * 3 + 3])

                        # Aplicar transformações e projeção
                        p1_2d, z1 = GL.project_vertex(p1)
                        p2_2d, z2 = GL.project_vertex(p2)
                        p3_2d, z3 = GL.project_vertex(p3)

                        colors_for_interpol = [emissiveColor, emissiveColor, emissiveColor]

                        if colorPerVertex and colorIndex:
                            c1 = np.array(color[color_indices[0] * 3: color_indices[0] * 3 + 3]) * 255
                            c2 = np.array(color[color_indices[j] * 3: color_indices[j] * 3 + 3]) * 255
                            c3 = np.array(color[color_indices[j + 1] * 3: color_indices[j + 1] * 3 + 3]) * 255

                            colors_for_interpol = [c1, c2, c3]

                        # Desenhar o triângulo
                        GL.rasterize_triangle([p1_2d, p2_2d, p3_2d], colors_for_interpol, z_coords=[z1, z2, z3])

                # Limpar a lista de índices para a próxima face
                face_indices = []
                color_indices = []
            else:
                face_indices.append(coordIndex[i])
                if colorIndex:
                    color_indices.append(colorIndex[i])
            i += 1


    @staticmethod
    def box(size, colors):
        """Função usada para renderizar Boxes."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry3D.html#Box
        # A função box é usada para desenhar paralelepípedos na cena. O Box é centrada no
        # (0, 0, 0) no sistema de coordenadas local e alinhado com os eixos de coordenadas
        # locais. O argumento size especifica as extensões da caixa ao longo dos eixos X, Y
        # e Z, respectivamente, e cada valor do tamanho deve ser maior que zero. Para desenha
        # essa caixa você vai provavelmente querer tesselar ela em triângulos, para isso
        # encontre os vértices e defina os triângulos.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("Box : size = {0}".format(size)) # imprime no terminal pontos
        print("Box : colors = {0}".format(colors)) # imprime no terminal as cores

        # Exemplo de desenho de um pixel branco na coordenada 10, 10
        gpu.GPU.draw_pixel([10, 10], gpu.GPU.RGB8, [255, 255, 255])  # altera pixel

    @staticmethod
    def sphere(radius, colors):
        """Função usada para renderizar Esferas."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry3D.html#Sphere
        # A função sphere é usada para desenhar esferas na cena. O esfera é centrada no
        # (0, 0, 0) no sistema de coordenadas local. O argumento radius especifica o
        # raio da esfera que está sendo criada. Para desenha essa esfera você vai
        # precisar tesselar ela em triângulos, para isso encontre os vértices e defina
        # os triângulos.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("Sphere : radius = {0}".format(radius)) # imprime no terminal o raio da esfera
        print("Sphere : colors = {0}".format(colors)) # imprime no terminal as cores

    @staticmethod
    def cone(bottomRadius, height, colors):
        """Função usada para renderizar Cones."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry3D.html#Cone
        # A função cone é usada para desenhar cones na cena. O cone é centrado no
        # (0, 0, 0) no sistema de coordenadas local. O argumento bottomRadius especifica o
        # raio da base do cone e o argumento height especifica a altura do cone.
        # O cone é alinhado com o eixo Y local. O cone é fechado por padrão na base.
        # Para desenha esse cone você vai precisar tesselar ele em triângulos, para isso
        # encontre os vértices e defina os triângulos.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("Cone : bottomRadius = {0}".format(bottomRadius)) # imprime no terminal o raio da base do cone
        print("Cone : height = {0}".format(height)) # imprime no terminal a altura do cone
        print("Cone : colors = {0}".format(colors)) # imprime no terminal as cores

    @staticmethod
    def cylinder(radius, height, colors):
        """Função usada para renderizar Cilindros."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry3D.html#Cylinder
        # A função cylinder é usada para desenhar cilindros na cena. O cilindro é centrado no
        # (0, 0, 0) no sistema de coordenadas local. O argumento radius especifica o
        # raio da base do cilindro e o argumento height especifica a altura do cilindro.
        # O cilindro é alinhado com o eixo Y local. O cilindro é fechado por padrão em ambas as extremidades.
        # Para desenha esse cilindro você vai precisar tesselar ele em triângulos, para isso
        # encontre os vértices e defina os triângulos.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("Cylinder : radius = {0}".format(radius)) # imprime no terminal o raio do cilindro
        print("Cylinder : height = {0}".format(height)) # imprime no terminal a altura do cilindro
        print("Cylinder : colors = {0}".format(colors)) # imprime no terminal as cores

    @staticmethod
    def navigationInfo(headlight):
        """Características físicas do avatar do visualizador e do modelo de visualização."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/navigation.html#NavigationInfo
        # O campo do headlight especifica se um navegador deve acender um luz direcional que
        # sempre aponta na direção que o usuário está olhando. Definir este campo como TRUE
        # faz com que o visualizador forneça sempre uma luz do ponto de vista do usuário.
        # A luz headlight deve ser direcional, ter intensidade = 1, cor = (1 1 1),
        # ambientIntensity = 0,0 e direção = (0 0 −1).

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("NavigationInfo : headlight = {0}".format(headlight)) # imprime no terminal

    @staticmethod
    def directionalLight(ambientIntensity, color, intensity, direction):
        """Luz direcional ou paralela."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/lighting.html#DirectionalLight
        # Define uma fonte de luz direcional que ilumina ao longo de raios paralelos
        # em um determinado vetor tridimensional. Possui os campos básicos ambientIntensity,
        # cor, intensidade. O campo de direção especifica o vetor de direção da iluminação
        # que emana da fonte de luz no sistema de coordenadas local. A luz é emitida ao
        # longo de raios paralelos de uma distância infinita.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("DirectionalLight : ambientIntensity = {0}".format(ambientIntensity))
        print("DirectionalLight : color = {0}".format(color)) # imprime no terminal
        print("DirectionalLight : intensity = {0}".format(intensity)) # imprime no terminal
        print("DirectionalLight : direction = {0}".format(direction)) # imprime no terminal

    @staticmethod
    def pointLight(ambientIntensity, color, intensity, location):
        """Luz pontual."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/lighting.html#PointLight
        # Fonte de luz pontual em um local 3D no sistema de coordenadas local. Uma fonte
        # de luz pontual emite luz igualmente em todas as direções; ou seja, é omnidirecional.
        # Possui os campos básicos ambientIntensity, cor, intensidade. Um nó PointLight ilumina
        # a geometria em um raio de sua localização. O campo do raio deve ser maior ou igual a
        # zero. A iluminação do nó PointLight diminui com a distância especificada.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("PointLight : ambientIntensity = {0}".format(ambientIntensity))
        print("PointLight : color = {0}".format(color)) # imprime no terminal
        print("PointLight : intensity = {0}".format(intensity)) # imprime no terminal
        print("PointLight : location = {0}".format(location)) # imprime no terminal

    @staticmethod
    def fog(visibilityRange, color):
        """Névoa."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/environmentalEffects.html#Fog
        # O nó Fog fornece uma maneira de simular efeitos atmosféricos combinando objetos
        # com a cor especificada pelo campo de cores com base nas distâncias dos
        # vários objetos ao visualizador. A visibilidadeRange especifica a distância no
        # sistema de coordenadas local na qual os objetos são totalmente obscurecidos
        # pela névoa. Os objetos localizados fora de visibilityRange do visualizador são
        # desenhados com uma cor de cor constante. Objetos muito próximos do visualizador
        # são muito pouco misturados com a cor do nevoeiro.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("Fog : color = {0}".format(color)) # imprime no terminal
        print("Fog : visibilityRange = {0}".format(visibilityRange))

    @staticmethod
    def timeSensor(cycleInterval, loop):
        """Gera eventos conforme o tempo passa."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/time.html#TimeSensor
        # Os nós TimeSensor podem ser usados para muitas finalidades, incluindo:
        # Condução de simulações e animações contínuas; Controlar atividades periódicas;
        # iniciar eventos de ocorrência única, como um despertador;
        # Se, no final de um ciclo, o valor do loop for FALSE, a execução é encerrada.
        # Por outro lado, se o loop for TRUE no final de um ciclo, um nó dependente do
        # tempo continua a execução no próximo ciclo. O ciclo de um nó TimeSensor dura
        # cycleInterval segundos. O valor de cycleInterval deve ser maior que zero.

        # Deve retornar a fração de tempo passada em fraction_changed

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("TimeSensor : cycleInterval = {0}".format(cycleInterval)) # imprime no terminal
        print("TimeSensor : loop = {0}".format(loop))

        # Esse método já está implementado para os alunos como exemplo
        epoch = time.time()  # time in seconds since the epoch as a floating point number.
        fraction_changed = (epoch % cycleInterval) / cycleInterval

        return fraction_changed

    @staticmethod
    def splinePositionInterpolator(set_fraction, key, keyValue, closed):
        """Interpola não linearmente entre uma lista de vetores 3D."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/interpolators.html#SplinePositionInterpolator
        # Interpola não linearmente entre uma lista de vetores 3D. O campo keyValue possui
        # uma lista com os valores a serem interpolados, key possui uma lista respectiva de chaves
        # dos valores em keyValue, a fração a ser interpolada vem de set_fraction que varia de
        # zeroa a um. O campo keyValue deve conter exatamente tantos vetores 3D quanto os
        # quadros-chave no key. O campo closed especifica se o interpolador deve tratar a malha
        # como fechada, com uma transições da última chave para a primeira chave. Se os keyValues
        # na primeira e na última chave não forem idênticos, o campo closed será ignorado.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("SplinePositionInterpolator : set_fraction = {0}".format(set_fraction))
        print("SplinePositionInterpolator : key = {0}".format(key)) # imprime no terminal
        print("SplinePositionInterpolator : keyValue = {0}".format(keyValue))
        print("SplinePositionInterpolator : closed = {0}".format(closed))

        # Abaixo está só um exemplo de como os dados podem ser calculados e transferidos
        value_changed = [0.0, 0.0, 0.0]
        
        return value_changed

    @staticmethod
    def orientationInterpolator(set_fraction, key, keyValue):
        """Interpola entre uma lista de valores de rotação especificos."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/interpolators.html#OrientationInterpolator
        # Interpola rotações são absolutas no espaço do objeto e, portanto, não são cumulativas.
        # Uma orientação representa a posição final de um objeto após a aplicação de uma rotação.
        # Um OrientationInterpolator interpola entre duas orientações calculando o caminho mais
        # curto na esfera unitária entre as duas orientações. A interpolação é linear em
        # comprimento de arco ao longo deste caminho. Os resultados são indefinidos se as duas
        # orientações forem diagonalmente opostas. O campo keyValue possui uma lista com os
        # valores a serem interpolados, key possui uma lista respectiva de chaves
        # dos valores em keyValue, a fração a ser interpolada vem de set_fraction que varia de
        # zeroa a um. O campo keyValue deve conter exatamente tantas rotações 3D quanto os
        # quadros-chave no key.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("OrientationInterpolator : set_fraction = {0}".format(set_fraction))
        print("OrientationInterpolator : key = {0}".format(key)) # imprime no terminal
        print("OrientationInterpolator : keyValue = {0}".format(keyValue))

        # Abaixo está só um exemplo de como os dados podem ser calculados e transferidos
        value_changed = [0, 0, 1, 0]

        return value_changed

    # Para o futuro (Não para versão atual do projeto.)
    def vertex_shader(self, shader):
        """Para no futuro implementar um vertex shader."""

    def fragment_shader(self, shader):
        """Para no futuro implementar um fragment shader."""
