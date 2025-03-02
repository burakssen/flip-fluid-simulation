#version 100
precision mediump float;

// Input vertex attributes (from vertex shader)
varying vec3 fragColor;
varying float fragDrawDisk;

// Output fragment color
// Note: In WebGL/GLSL 100, we don't declare outputs with 'out'

void main()
{
    if (fragDrawDisk > 0.9) {  // Using approx comparison for floating point
        float rx = 0.5 - gl_PointCoord.x;
        float ry = 0.5 - gl_PointCoord.y;
        float r2 = rx * rx + ry * ry;
        
        if (r2 > 0.25)
            discard;
    }
    
    gl_FragColor = vec4(0.0, 1.0, 1.0, 1.0);
}