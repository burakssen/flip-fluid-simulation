#include <raylib.h>
#include <rlgl.h>
#include <raymath.h>
#include "FlipFluid.h"
#include <vector>
#include <string>
#include <memory>
#include <cmath>

// Scene configuration
struct FluidScene
{
    // Simulation parameters
    float gravity = -9.81f;
    float dt = 1.0f / 60.0f;
    float flipRatio = 0.9f;
    int numPressureIters = 100;
    int numParticleIters = 2;
    float overRelaxation = 1.9f;

    // State tracking
    int frameNr = 0;
    bool paused = true;
    bool showParticles = true;

    // Obstacle properties
    float obstacleX = 0.0f;
    float obstacleY = 0.0f;
    float obstacleRadius = 0.15f;
    float obstacleVelX = 0.0f;
    float obstacleVelY = 0.0f;
    bool showObstacle = true;

    // Fluid simulation
    std::unique_ptr<FlipFluid> fluid;
};

// Global variables
FluidScene scene;
float screenWidth = 1280;
float screenHeight = 720;
float cScale;
bool mouseDown = false;
Vector2 lastMousePos = {0};
Camera2D camera = {0};

// Shader variables
Shader particleShader;
int particlePointSizeLoc;
int domainSizeLoc;
int drawDiskLoc;

// Forward declarations
void InitializeScene(float width, float height);
void SetObstacle(float x, float y, bool reset);
void DrawScene();
void SimulateStep();
void HandleInput();
void RenderFrame();

int main()
{
    InitWindow(screenWidth, screenHeight, "FLIP Fluid Simulation");

    // Setup camera
    camera.offset = {screenWidth / 2.0f, screenHeight / 2.0f};
    camera.target = {screenWidth / 2.0f, screenHeight / 2.0f};
    camera.zoom = 1.0f;

    // Calculate simulation dimensions
    float simHeight = 3.0f;
    cScale = screenHeight / simHeight;
    float simWidth = screenWidth / cScale;

    InitializeScene(simWidth, simHeight);

    // Load shaders
    particleShader = LoadShader("../shaders/point.vs", "../shaders/point.fs");
    particlePointSizeLoc = GetShaderLocation(particleShader, "pointSize");
    domainSizeLoc = GetShaderLocation(particleShader, "domainSize");
    drawDiskLoc = GetShaderLocation(particleShader, "drawDisk");

    // Main game loop
    while (!WindowShouldClose())
    {
        HandleInput();
        if (!scene.paused)
            SimulateStep();
        RenderFrame();
    }

    UnloadShader(particleShader);
    CloseWindow();
    return 0;
}

void HandleInput()
{
    // Toggle pause state
    if (IsKeyPressed(KEY_P))
        scene.paused = !scene.paused;

    // Single step simulation
    if (IsKeyPressed(KEY_M))
    {
        scene.paused = false;
        SimulateStep();
        scene.paused = true;
    }

    // Mouse interaction for obstacle control
    if (IsMouseButtonDown(MOUSE_LEFT_BUTTON))
    {
        Vector2 mousePos = GetMousePosition();
        float x = mousePos.x / cScale;
        float y = (screenHeight - mousePos.y) / cScale;

        if (!mouseDown)
        {
            SetObstacle(x, y, true);
            mouseDown = true;
            scene.paused = false;
        }
        else
        {
            SetObstacle(x, y, false);
        }

        lastMousePos = mousePos;
    }

    if (IsMouseButtonReleased(MOUSE_LEFT_BUTTON))
    {
        mouseDown = false;
        scene.obstacleVelX = 0.0f;
        scene.obstacleVelY = 0.0f;
    }
}

void RenderFrame()
{
    BeginDrawing();
    ClearBackground(BLACK);

    BeginMode2D(camera);
    DrawScene();
    EndMode2D();

    // UI text
    DrawText(scene.paused ? "PAUSED (Press P to resume)" : "Running (Press P to pause)", 10, 10, 20, WHITE);
    DrawFPS(screenWidth - 100, 10);

    EndDrawing();
}

void InitializeScene(float width, float height)
{
    // Tank dimensions and grid resolution
    float tankHeight = height;
    float tankWidth = width;
    float h = tankHeight / 100.0f;
    float density = 1000.0f;

    // Water volume parameters (relative to tank size)
    float relWaterHeight = 0.8f;
    float relWaterWidth = 0.8f;

    // Particle parameters and spacing
    float r = 0.3f * h;
    float dx = 2.0f * r;
    float dy = (sqrtf(3.0f) / 2.0f) * dx;

    // Calculate max particle count
    int numX = static_cast<int>(std::floor((relWaterWidth * tankWidth - 2.0f * h - 2.0f * r) / dx));
    int numY = static_cast<int>(std::floor((relWaterHeight * tankHeight - 2.0f * h - 2.0f * r) / dy));
    int maxParticles = numX * numY;

    // Create fluid simulator
    scene.fluid = std::make_unique<FlipFluid>(density, tankWidth, tankHeight, h, r, maxParticles);

    // Initialize dam break scenario
    scene.fluid->numParticles = numX * numY;
    int p = 0;
    for (int i = 0; i < numX; i++)
    {
        for (int j = 0; j < numY; j++)
        {
            scene.fluid->particlePos[p++] = h + r + dx * i + ((j % 2 == 0) ? 0.0f : r);
            scene.fluid->particlePos[p++] = h + r + dy * j;
        }
    }

    // Setup boundaries (solid cells)
    int n = scene.fluid->fNumY;
    for (int i = 0; i < scene.fluid->fNumX; i++)
    {
        for (int j = 0; j < scene.fluid->fNumY; j++)
        {
            float s = (i == 0 || i == scene.fluid->fNumX - 1 || j == 0) ? 0.0f : 1.0f;
            scene.fluid->s[i * n + j] = s;
        }
    }

    // Set the initial obstacle position (center of the tank)
    SetObstacle(width / 2.0f, height / 2.0f, true);
}

void SetObstacle(float x, float y, bool reset)
{
    if (!reset)
    {
        scene.obstacleVelX = (x - scene.obstacleX) / scene.dt;
        scene.obstacleVelY = (y - scene.obstacleY) / scene.dt;
    }
    scene.obstacleX = x;
    scene.obstacleY = y;
    scene.showObstacle = true;
}

void SimulateStep()
{
    scene.fluid->simulate(
        scene.dt,
        scene.gravity,
        scene.flipRatio,
        scene.numPressureIters,
        scene.numParticleIters,
        scene.overRelaxation,
        scene.obstacleX,
        scene.obstacleY,
        scene.obstacleRadius,
        scene.obstacleVelX,
        scene.obstacleVelY);
    scene.frameNr++;
}

void DrawScene()
{
    // Render particles using the shader and RL_LINES
    if (scene.showParticles)
    {
        BeginShaderMode(particleShader);

        // Pass uniforms to shader
        float pointSize = ((2.0f * scene.fluid->particleRadius) / (screenWidth / cScale)) * screenWidth;
        SetShaderValue(particleShader, particlePointSizeLoc, &pointSize, SHADER_UNIFORM_FLOAT);

        float domainSize[2] = {screenWidth / cScale, screenHeight / cScale};
        SetShaderValue(particleShader, domainSizeLoc, domainSize, SHADER_UNIFORM_VEC2);

        int drawDisk = scene.showObstacle ? 1 : 0;
        SetShaderValue(particleShader, drawDiskLoc, &drawDisk, SHADER_UNIFORM_INT);

        const auto &particlePos = scene.fluid->particlePos;
        const auto &particleColor = scene.fluid->particleColor;
        int numParticles = scene.fluid->numParticles;

        // Draw particles as very short lines (point simulation)
        rlBegin(RL_LINES);
        for (int i = 0; i < numParticles; i++)
        {
            float x = particlePos[2 * i];
            float y = particlePos[2 * i + 1];
            float screenX = x * cScale;
            float screenY = screenHeight - y * cScale;

            Color color = {
                static_cast<unsigned char>(particleColor[3 * i] * 255),
                static_cast<unsigned char>(particleColor[3 * i + 1] * 255),
                static_cast<unsigned char>(particleColor[3 * i + 2] * 255),
                255};

            rlColor4ub(color.r, color.g, color.b, color.a);
            rlVertex3f(screenX, screenY, 0.0f);
            rlVertex3f(screenX + 1.0f, screenY, 0.0f);
        }
        rlEnd();

        EndShaderMode();
    }

    // Draw the obstacle
    if (scene.showObstacle)
    {
        float screenX = scene.obstacleX * cScale;
        float screenY = screenHeight - scene.obstacleY * cScale;
        float screenRadius = scene.obstacleRadius * cScale;
        DrawCircle(static_cast<int>(screenX), static_cast<int>(screenY),
                   static_cast<int>(screenRadius), RED);
    }
}