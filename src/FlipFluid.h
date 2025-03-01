#pragma once

#include <vector>
#include <algorithm>
#include <cmath>

#include <raymath.h>

enum CellType
{
    FLUID_CELL = 0,
    AIR_CELL = 1,
    SOLID_CELL = 2
};

class FlipFluid
{
public:
    FlipFluid(float density, float width, float height, float spacing, float particleRadius, int maxParticles);

private:
    void integrateParticles(float dt, float gravity);
    void pushParticlesApart(int numIters);
    void handleParticleCollisions(float obstacleX, float obstacleY, float obstacleRadius, float obstacleVelX, float obstacleVelY);
    void updateParticleDensity();
    void transferVelocities(bool toGrid, float flipRatio);
    void solveIncompressibility(int numIters, float dt, float overRelaxation);
    void updateParticleColors();
    void setSciColor(int cellNr, float val, float minVal, float maxVal);
    void updateCellColors();

public:
    void simulate(float dt, float gravity, float flipRatio, int numPressureIters, int numParticleIters, float overRelaxation, float obstacleX, float obstacleY, float obstacleRadius, float obstacleVelX, float obstacleVelY)
    {
        const int numSubSteps = 1;
        const float sdt = dt / numSubSteps;
        for (int step = 0; step < numSubSteps; step++)
        {
            this->integrateParticles(sdt, gravity);
            this->pushParticlesApart(numParticleIters);
            this->handleParticleCollisions(obstacleX, obstacleY, obstacleRadius, obstacleVelX, obstacleVelY);
            this->transferVelocities(true, 0.0f);
            this->updateParticleDensity();
            this->solveIncompressibility(numPressureIters, sdt, overRelaxation);
            this->transferVelocities(false, flipRatio);
        }

        this->updateParticleColors();
        this->updateCellColors();
    }

public:
    float density;
    Vector2 fNum;
    int fNumX;
    int fNumY;
    float h;
    float fInvSpacing;
    int fNumCells;
    std::vector<float> u;
    std::vector<float> v;
    std::vector<float> du;
    std::vector<float> dv;
    std::vector<float> prevU;
    std::vector<float> prevV;
    std::vector<float> p;
    std::vector<float> s;
    std::vector<CellType> cellType;
    std::vector<float> cellColor;

    int maxParticles;
    std::vector<float> particlePos;
    std::vector<float> particleColor;
    std::vector<float> particleVel;
    std::vector<float> particleDensity;
    float particleRestDensity;
    float particleRadius;
    float pInvSpacing;
    int pNumX;
    int pNumY;
    int pNumCells;

    std::vector<int> numCellParticles;
    std::vector<int> firstCellParticle;
    std::vector<int> cellParticleIds;

    int numParticles;
};