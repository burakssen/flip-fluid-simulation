#include "FlipFluid.h"

FlipFluid::FlipFluid(float density, float width, float height, float spacing, float particleRadius, int maxParticles)
{
    this->density = density;
    this->fNumX = std::floor(width / spacing) + 1;
    this->fNumY = std::floor(height / spacing) + 1;
    this->h = std::max(width / this->fNumX, height / this->fNumY);
    this->fInvSpacing = 1.0f / this->h;
    this->fNumCells = this->fNumX * this->fNumY;

    this->u.resize(this->fNumCells);
    this->v.resize(this->fNumCells);
    this->du.resize(this->fNumCells);
    this->dv.resize(this->fNumCells);
    this->prevU.resize(this->fNumCells);
    this->prevV.resize(this->fNumCells);
    this->p.resize(this->fNumCells);
    this->s.resize(this->fNumCells);
    this->cellType.resize(this->fNumCells);
    this->cellColor.resize(3 * this->fNumCells);

    this->maxParticles = maxParticles;
    this->particlePos.resize(2 * this->maxParticles);
    this->particleColor.resize(3 * this->maxParticles);

    for (int i = 0; i < this->maxParticles; i++)
        this->particleColor[3 * i + 2.0] = 1.0f;

    this->particleVel.resize(2 * this->maxParticles);
    this->particleDensity.resize(this->fNumCells);
    this->particleRestDensity = 0.0f;

    this->particleRadius = particleRadius;
    this->pInvSpacing = 1.0f / (2.2f * particleRadius);
    this->pNumX = std::floor(width * this->pInvSpacing) + 1;
    this->pNumY = std::floor(height * this->pInvSpacing) + 1;
    this->pNumCells = this->pNumX * this->pNumY;

    this->numCellParticles.resize(this->pNumCells);
    this->firstCellParticle.resize(this->pNumCells + 1);
    this->cellParticleIds.resize(this->maxParticles);

    this->numParticles = 0;
}

void FlipFluid::integrateParticles(float dt, float gravity)
{
    for (int i = 0; i < this->numParticles; i++)
    {
        float vx = this->particleVel[2 * i];
        float vy = this->particleVel[2 * i + 1] + dt * gravity;

        this->particleVel[2 * i + 1] = vy;
        this->particlePos[2 * i] += vx * dt;
        this->particlePos[2 * i + 1] += vy * dt;
    }
}

void FlipFluid::pushParticlesApart(int numIters)
{
    this->numCellParticles.assign(this->pNumCells, 0);

    for (int i = 0; i < this->numParticles; i++)
    {
        const float x = this->particlePos[2 * i];
        const float y = this->particlePos[2 * i + 1];
        const int xi = std::clamp((int)std::floor(x * this->pInvSpacing), 0, this->pNumX - 1);
        const int yi = std::clamp((int)std::floor(y * this->pInvSpacing), 0, this->pNumY - 1);
        const int cellNr = xi * this->pNumY + yi;
        this->numCellParticles[cellNr]++;
    }

    int sum = 0;
    for (int i = 0; i < this->pNumCells; i++)
    {
        sum += this->numCellParticles[i];
        this->firstCellParticle[i] = sum;
    }

    this->firstCellParticle[this->pNumCells] = sum;

    for (int i = 0; i < this->numParticles; i++)
    {
        const float x = this->particlePos[2 * i];
        const float y = this->particlePos[2 * i + 1];
        const int xi = std::clamp((int)std::floor(x * this->pInvSpacing), 0, this->pNumX - 1);
        const int yi = std::clamp((int)std::floor(y * this->pInvSpacing), 0, this->pNumY - 1);
        const int cellNr = xi * this->pNumY + yi;
        this->firstCellParticle[cellNr]--;
        this->cellParticleIds[this->firstCellParticle[cellNr]] = i;
    }

    float minDist = 2.0f * this->particleRadius;
    float minDist2 = minDist * minDist;

    for (int iter = 0; iter < numIters; iter++)
    {
        for (int i = 0; i < this->numParticles; i++)
        {
            const float px = this->particlePos[2 * i];
            const float py = this->particlePos[2 * i + 1];

            const int pxi = std::floor(px * this->pInvSpacing);
            const int pyi = std::floor(py * this->pInvSpacing);
            const int x0 = std::max(pxi - 1, 0);
            const int y0 = std::max(pyi - 1, 0);
            const int x1 = std::min(pxi + 1, this->pNumX - 1);
            const int y1 = std::min(pyi + 1, this->pNumY - 1);

            for (int xi = x0; xi <= x1; xi++)
            {
                for (int yi = y0; yi <= y1; yi++)
                {
                    const int cellNr = xi * this->pNumY + yi;
                    const int start = this->firstCellParticle[cellNr];
                    const int end = this->firstCellParticle[cellNr + 1];

                    for (int j = start; j < end; j++)
                    {
                        const int id = this->cellParticleIds[j];
                        if (id == i)
                            continue;
                        const float qx = this->particlePos[2 * id];
                        const float qy = this->particlePos[2 * id + 1];

                        const float dx = qx - px;
                        const float dy = qy - py;
                        const float d2 = dx * dx + dy * dy;
                        if (d2 > minDist2 || d2 == 0.0f)
                            continue;

                        const float d = std::sqrt(d2);
                        const float s = (0.5f * (minDist - d)) / d;
                        const float offsetX = dx * s;
                        const float offsetY = dy * s;

                        this->particlePos[2 * i] -= offsetX;
                        this->particlePos[2 * i + 1] -= offsetY;
                        this->particlePos[2 * id] += offsetX;
                        this->particlePos[2 * id + 1] += offsetY;
                    }
                }
            }
        }
    }
}

void FlipFluid::handleParticleCollisions(float obstacleX, float obstacleY, float obstacleRadius, float obstacleVelX, float obstacleVelY)
{
    const float h = 1.0f / this->fInvSpacing;
    const float minX = h + this->particleRadius;
    const float maxX = (this->fNumX - 1) * h - this->particleRadius;
    const float minY = minX;
    const float maxY = maxX;

    const float minDist2 = (obstacleRadius + this->particleRadius) * (obstacleRadius + this->particleRadius);

    for (int i = 0; i < this->numParticles; i++)
    {
        float x = this->particlePos[2 * i];
        float y = this->particlePos[2 * i + 1];

        const float dx = x - obstacleX;
        const float dy = y - obstacleY;
        const float d2 = dx * dx + dy * dy;

        if (d2 < minDist2)
        {
            this->particleVel[2 * i] = obstacleVelX;
            this->particleVel[2 * i + 1] = obstacleVelY;
        }

        if (x < minX || x > maxX)
        {
            this->particleVel[2 * i] = 0.0f;
            x = std::max(minX, std::min(x, maxX));
        }

        if (y < minY || y > maxY)
        {
            this->particleVel[2 * i + 1] = 0.0f;
            y = std::max(minY, std::min(y, maxY));
        }

        this->particlePos[2 * i] = x;
        this->particlePos[2 * i + 1] = y;
    }
}

void FlipFluid::updateParticleDensity()
{
    const float n = this->fNumY;
    const float h1 = this->fInvSpacing;
    const float h2 = 0.5f * h1;

    this->particleDensity.assign(this->fNumCells, 0.0f);
    for (int i = 0; i < this->numParticles; i++)
    {
        float x = std::clamp(this->particlePos[2 * i], this->h, (this->fNumX - 1) * this->h);
        float y = std::clamp(this->particlePos[2 * i + 1], this->h, (this->fNumY - 1) * this->h);

        const int x0 = std::floor((x - h2) * h1);
        const int y0 = std::floor((y - h2) * h1);
        const int x1 = std::min(x0 + 1, this->fNumX - 2);
        const int y1 = std::min(y0 + 1, this->fNumY - 2);

        const float tx = (x - h2 - x0 * this->h) * h1;
        const float ty = (y - h2 - y0 * this->h) * h1;
        const float sx = 1.0f - tx;
        const float sy = 1.0f - ty;

        this->particleDensity[x0 * n + y0] += sx * sy;
        this->particleDensity[x1 * n + y0] += tx * sy;
        this->particleDensity[x1 * n + y1] += tx * ty;
        this->particleDensity[x0 * n + y1] += sx * ty;
    }

    if (this->particleRestDensity == 0.0f)
    {
        float sum = 0.0f;
        int numFluidCells = 0;
        for (int i = 0; i < this->fNumCells; i++)
        {
            if (this->cellType[i] == CellType::FLUID_CELL)
            {
                sum += this->particleDensity[i];
                numFluidCells++;
            }
        }

        if (numFluidCells > 0)
        {
            this->particleRestDensity = sum / numFluidCells;
        }
    }
}

void FlipFluid::transferVelocities(bool toGrid, float flipRatio)
{
    const int n = this->fNumY;
    const float h = this->h;
    const float h1 = this->fInvSpacing;
    const float h2 = 0.5f * h;

    if (toGrid)
    {
        this->prevU = this->u;
        this->prevV = this->v;
        this->du.assign(this->fNumCells, 0.0f);
        this->dv.assign(this->fNumCells, 0.0f);
        this->u.assign(this->fNumCells, 0.0f);
        this->v.assign(this->fNumCells, 0.0f);

        for (int i = 0; i < this->fNumCells; i++)
        {
            this->cellType[i] = this->s[i] == 0.0f ? CellType::SOLID_CELL : CellType::AIR_CELL;
        }

        for (int i = 0; i < this->numParticles; i++)
        {
            const float x = this->particlePos[2 * i];
            const float y = this->particlePos[2 * i + 1];
            const int xi = std::clamp((int)std::floor(x * h1), 0, this->fNumX - 1);
            const int yi = std::clamp((int)std::floor(y * h1), 0, this->fNumY - 1);
            const int cellNr = xi * n + yi;

            if (this->cellType[cellNr] == CellType::AIR_CELL)
            {
                this->cellType[cellNr] = CellType::FLUID_CELL;
            }
        }
    }

    for (int component = 0; component < 2; component++)
    {
        const float dx = component == 0 ? 0.0f : h2;
        const float dy = component == 0 ? h2 : 0.0f;
        std::vector<float> &f = component == 0 ? this->u : this->v;
        std::vector<float> &prevF = component == 0 ? this->prevU : this->prevV;
        std::vector<float> &d = component == 0 ? this->du : this->dv;
        int offset = component == 0 ? n : 1;

        for (int i = 0; i < this->numParticles; i++)
        {
            float x = this->particlePos[2 * i];
            float y = this->particlePos[2 * i + 1];

            x = std::clamp(x, h, (this->fNumX - 1) * h);
            y = std::clamp(y, h, (this->fNumY - 1) * h);

            const int x0 = std::min((int)std::floor((x - dx) * h1), this->fNumX - 2);
            const float tx = (x - dx - x0 * h) * h1;
            const int x1 = std::min(x0 + 1, this->fNumX - 1);

            const int y0 = std::min((int)std::floor((y - dy) * h1), this->fNumY - 2);
            const float ty = (y - dy - y0 * h) * h1;
            const int y1 = std::min(y0 + 1, this->fNumY - 2);

            const float sx = 1.0f - tx;
            const float sy = 1.0f - ty;

            const float d0 = sx * sy;
            const float d1 = tx * sy;
            const float d2 = tx * ty;
            const float d3 = sx * ty;

            const int nr0 = x0 * n + y0;
            const int nr1 = x1 * n + y0;
            const int nr2 = x1 * n + y1;
            const int nr3 = x0 * n + y1;

            if (toGrid)
            {
                const float pv = this->particleVel[2 * i + component];

                f[nr0] += d0 * pv;
                f[nr1] += d1 * pv;
                f[nr2] += d2 * pv;
                f[nr3] += d3 * pv;
                d[nr0] += d0;
                d[nr1] += d1;
                d[nr2] += d2;
                d[nr3] += d3;
            }
            else
            {
                const float valid0 = this->cellType[nr0] != CellType::AIR_CELL || this->cellType[nr0 - offset] != CellType::AIR_CELL ? 1.0f : 0.0f;
                const float valid1 = this->cellType[nr1] != CellType::AIR_CELL || this->cellType[nr1 - offset] != CellType::AIR_CELL ? 1.0f : 0.0f;
                const float valid2 = this->cellType[nr2] != CellType::AIR_CELL || this->cellType[nr2 - offset] != CellType::AIR_CELL ? 1.0f : 0.0f;
                const float valid3 = this->cellType[nr3] != CellType::AIR_CELL || this->cellType[nr3 - offset] != CellType::AIR_CELL ? 1.0f : 0.0f;

                const float totalWeight = valid0 * d0 + valid1 * d1 + valid2 * d2 + valid3 * d3;

                if (totalWeight > 0.0f)
                {
                    const float v = this->particleVel[2 * i + component];
                    const float picV = (valid0 * d0 * f[nr0] + valid1 * d1 * f[nr1] + valid2 * d2 * f[nr2] + valid3 * d3 * f[nr3]) / totalWeight;
                    const float corr = (valid0 * d0 * (f[nr0] - prevF[nr0]) + valid1 * d1 * (f[nr1] - prevF[nr1]) + valid2 * d2 * (f[nr2] - prevF[nr2]) + valid3 * d3 * (f[nr3] - prevF[nr3])) / totalWeight;

                    this->particleVel[2 * i + component] = (1.0f - flipRatio) * picV + flipRatio * (v + corr);
                }
            }
        }

        if (toGrid)
        {
            for (int i = 0; i < f.size(); i++)
            {
                if (d[i] > 0.0f)
                {
                    f[i] /= d[i];
                }
            }

            for (int i = 0; i < this->fNumX; i++)
            {
                for (int j = 0; j < this->fNumY; j++)
                {
                    const int idx = i * n + j;
                    const bool solid = this->cellType[idx] == CellType::SOLID_CELL;

                    if (solid || (i > 0 && this->cellType[(i - 1) * n + j] == CellType::SOLID_CELL))
                    {
                        this->u[idx] = this->prevU[idx];
                    }

                    if (solid || (j > 0 && this->cellType[i * n + j - 1] == CellType::SOLID_CELL))
                    {
                        this->v[idx] = this->prevV[idx];
                    }
                }
            }
        }
    }
}

void FlipFluid::solveIncompressibility(int numIters, float dt, float overRelaxation)
{
    this->p.assign(this->fNumCells, 0.0f);
    this->prevU = this->u;
    this->prevV = this->v;

    const int n = this->fNumY;
    const float cp = (this->density * this->h) / dt;

    for (int iter = 0; iter < numIters; iter++)
    {
        for (int i = 1; i < this->fNumX - 1; i++)
        {
            for (int j = 1; j < this->fNumY - 1; j++)
            {
                const int center = i * n + j;
                if (this->cellType[center] != CellType::FLUID_CELL)
                    continue;

                const int left = (i - 1) * n + j;
                const int right = (i + 1) * n + j;
                const int bottom = i * n + j - 1;
                const int top = i * n + j + 1;

                const float sx0 = this->s[left];
                const float sx1 = this->s[right];
                const float sy0 = this->s[bottom];
                const float sy1 = this->s[top];

                const float s = sx0 + sx1 + sy0 + sy1;
                if (s == 0.0f)
                    continue;

                float div = this->u[right] - this->u[center] + this->v[top] - this->v[center];

                if (this->particleRestDensity > 0.0f)
                {
                    const float compression = this->particleDensity[center] - this->particleRestDensity;
                    if (compression > 0.0f)
                    {
                        div -= compression;
                    }
                }

                const float p = (-div / s) * overRelaxation;

                this->p[center] += cp * p;
                this->u[center] -= sx0 * p;
                this->u[right] += sx1 * p;
                this->v[center] -= sy0 * p;
                this->v[top] += sy1 * p;
            }
        }
    }
}

void FlipFluid::updateParticleColors()
{
    const float h1 = this->fInvSpacing;
    const float s = 0.01f;
    const float d0 = this->particleRestDensity;
    const bool hasRestDensity = d0 > 0.0f;

    for (int i = 0; i < this->numParticles; i++)
    {
        this->particleColor[3 * i] = std::clamp(
            this->particleColor[3 * i] - s,
            0.0f,
            1.0f);
        this->particleColor[3 * i + 1] = std::clamp(
            this->particleColor[3 * i + 1] - s,
            0.0f,
            1.0f);
        this->particleColor[3 * i + 2] = std::clamp(
            this->particleColor[3 * i + 2] + s,
            0.0f,
            1.0f);

        if (hasRestDensity)
        {
            const float x = this->particlePos[2 * i];
            const float y = this->particlePos[2 * i + 1];
            const int xi = std::clamp((int)std::floor(x * h1), 1, this->fNumX - 1);
            const int yi = std::clamp((int)std::floor(y * h1), 1, this->fNumY - 1);
            const int cellNr = xi * this->fNumY + yi;
            const float relDensity = this->particleDensity[cellNr] / d0;

            if (relDensity < 0.7)
            {
                this->particleColor[3 * i] = 0.8f;
                this->particleColor[3 * i + 1] = 0.8f;
                this->particleColor[3 * i + 2] = 1.0f;
            }
        }
    }
}

void FlipFluid::setSciColor(int cellNr, float value, float minVal, float maxVal)
{
    float val = std::min(std::max(value, minVal), maxVal - 0.0001f);

    const float d = maxVal - minVal;
    val = d == 0.0f ? 0.5f : (val - minVal) / d;

    const float m = 0.25f;
    const int num = std::floor(val / m);
    const float s = (val - num * m) / m;

    float r, g, b;
    switch (num)
    {
    case 0:
        r = 0.0f;
        g = s;
        b = 1.0f;
        break;
    case 1:
        r = 0.0f;
        g = 1.0f;
        b = 1.0f - s;
        break;
    case 2:
        r = s;
        g = 1.0f;
        b = 0.0f;
        break;
    case 3:
        r = 1.0f;
        g = 1.0f - s;
        b = 0.0f;
        break;
    default:
        break;
    }

    const int colorIdx = 3 * cellNr;
    this->cellColor[colorIdx] = r;
    this->cellColor[colorIdx + 1] = g;
    this->cellColor[colorIdx + 2] = b;
}

void FlipFluid::updateCellColors()
{
    this->cellColor.assign(3 * this->fNumCells, 0.0f);

    const bool hasRestDensity = this->particleRestDensity > 0.0f;

    for (int i = 0; i < this->fNumCells; i++)
    {
        const CellType cellType = this->cellType[i];
        const int colorIdx = 3 * i;

        if (cellType == CellType::SOLID_CELL)
        {
            this->cellColor[colorIdx] = 0.5f;
            this->cellColor[colorIdx + 1] = 0.5f;
            this->cellColor[colorIdx + 2] = 0.5f;
        }
        else if (cellType == CellType::FLUID_CELL)
        {
            float d = this->particleDensity[i];

            if (hasRestDensity)
            {
                d /= this->particleRestDensity;
            }

            this->setSciColor(i, d, 0.0f, 2.0f);
        }
    }
}
