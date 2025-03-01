#version 330

// Input vertex attributes
in vec3 vertexPosition;
in vec2 vertexTexCoord;
in vec3 vertexColor;

// Output vertex attributes
out vec3 fragColor;
out float fragDrawDisk;

// Uniform locations
uniform mat4 mvp;
uniform vec2 domainSize;
uniform float pointSize;
uniform float drawDisk;

void main()
{
    // Position with custom transform (you may want to use mvp instead depending on your needs)
    vec4 position = mvp * vec4(vertexPosition.xy, 0.0, 1.0);
    
    // Alternative using raylib's mvp matrix:
    // vec4 position = mvp * vec4(vertexPosition, 1.0);
    
    gl_Position = position;
    gl_PointSize = pointSize;
    
    // Pass to fragment shader
    fragColor = vertexColor;
    fragDrawDisk = drawDisk;
}