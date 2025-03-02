#version 330

// Input vertex attributes (from vertex shader)
in vec3 fragColor;
in float fragDrawDisk;

// Output fragment color
out vec4 finalColor;

void main()
{
    if (fragDrawDisk == 1.0) {
        float rx = 0.5 - gl_PointCoord.x;
        float ry = 0.5 - gl_PointCoord.y;
        float r2 = rx * rx + ry * ry;
        if (r2 > 0.25)
            discard;
    }
    
    finalColor = vec4(0.0, 1.0, 1.0, 1.0);
}