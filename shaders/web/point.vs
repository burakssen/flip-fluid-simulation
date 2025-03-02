#version 100

// Input vertex attributes
attribute vec3 vertexPosition;
attribute vec2 vertexTexCoord;
attribute vec3 vertexColor;

// Output vertex attributes (to fragment shader)
varying vec3 fragColor;
varying float fragDrawDisk;

// Uniform locations
uniform mat4 mvp;
uniform vec2 domainSize;
uniform float pointSize;
uniform float drawDisk;

void main()
{
    // Position with custom transform
    vec4 position = mvp * vec4(vertexPosition.xy, 0.0, 1.0);
    
    // Alternative using raylib's mvp matrix:
    // vec4 position = mvp * vec4(vertexPosition, 1.0);
    
    gl_Position = position;
    gl_PointSize = pointSize;
    
    // Pass to fragment shader
    fragColor = vertexColor;
    fragDrawDisk = drawDisk;
}