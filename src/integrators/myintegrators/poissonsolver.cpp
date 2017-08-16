#include "poissonsolver.h"

#include <iostream>

namespace ps {

void ApplyATA(const double alpha,
              const Image3 &input,
              const Image3 &irlsWeight,
              Image3 &output) {
    const int pixelWidth = input.pixelWidth;
    const int pixelHeight = input.pixelHeight;

    for (int y = 0; y < pixelHeight; y++) {
        for (int x = 0; x < pixelWidth; x++) {
            output.At(x, y) = alpha * input.At(x, y) * irlsWeight.At(x, y)[0];

            Vector3 fxx = Vector3::Zero();
            if (x == 0) {
                fxx = (input.At(x + 1, y) - input.At(x, y)) * irlsWeight.At(x, y)[1];
            } else if (x < pixelWidth - 1) {
                fxx = (input.At(x + 1, y) - input.At(x, y)) * irlsWeight.At(x, y)[1] +
                    (input.At(x - 1, y) - input.At(x, y)) * irlsWeight.At(x - 1, y)[1];
            } else {
                fxx = (-input.At(x, y)) * irlsWeight.At(x, y)[1] +
                    (-input.At(x, y) + input.At(x - 1, y)) * irlsWeight.At(x - 1, y)[1];
            }
            output.At(x, y) -= fxx;
            
            Vector3 fyy = Vector3::Zero();
            if (y == 0) {
                fyy = (input.At(x, y + 1) - input.At(x, y)) * irlsWeight.At(x, y)[2];
            } else if (y < pixelHeight - 1) {
                fyy = (input.At(x, y + 1) - input.At(x, y)) * irlsWeight.At(x, y)[2] +
                    (input.At(x, y - 1) - input.At(x, y)) * irlsWeight.At(x, y - 1)[2];
            } else {
                fyy = (-input.At(x, y)) * irlsWeight.At(x, y)[2] +
                    (-input.At(x, y) + input.At(x, y - 1)) * irlsWeight.At(x, y - 1)[2];
            }
            output.At(x, y) -= fyy;
        }
    }
}

void ComputeATb(const double alpha,
                const Image3 &throughputImage,
                const Image3 &dxImage,
                const Image3 &dyImage,
                const Image3 &irlsWeight,
                Image3 &atb) {
    const int pixelWidth = throughputImage.pixelWidth;
    const int pixelHeight = throughputImage.pixelHeight;
    for (int y = 0; y < pixelHeight; y++) {
        for (int x = 0; x < pixelWidth; x++) {
            atb.At(x, y) = alpha * throughputImage.At(x, y) * irlsWeight.At(x, y)[0];

            Vector3 dxx = Vector3::Zero();
            if (x == 0) {
                dxx = dxImage.At(x, y) * irlsWeight.At(x, y)[1];
            } else {
                dxx = dxImage.At(x, y) * irlsWeight.At(x, y)[1] -
                    dxImage.At(x - 1, y) * irlsWeight.At(x - 1, y)[1];
            }
            atb.At(x, y) -= dxx;

            Vector3 dyy = Vector3::Zero();
            if (y == 0) {
                dyy = dyImage.At(x, y) * irlsWeight.At(x, y)[2];
            } else {
                dyy = dyImage.At(x, y) * irlsWeight.At(x, y)[2] -
                    dyImage.At(x, y - 1) * irlsWeight.At(x, y - 1)[2];
            }
            atb.At(x, y) -= dyy;
        }
    }
}


void Add(const Image3 &a, const Image3 &b, Image3 &out) {
    const int pixelWidth = a.pixelWidth;
    const int pixelHeight = a.pixelHeight;
    for (int y = 0; y < pixelHeight; y++) {
        for (int x = 0; x < pixelWidth; x++) {
            out.At(x, y) = a.At(x, y) + b.At(x, y);
        }
    }
}

void AddWithMultiplier(const Image3 &a, const Vector3 &multiplier, const Image3 &b, Image3 &out) {
    const int pixelWidth = a.pixelWidth;
    const int pixelHeight = a.pixelHeight;
    for (int y = 0; y < pixelHeight; y++) {
        for (int x = 0; x < pixelWidth; x++) {
            out.At(x, y) = a.At(x, y) + multiplier.cwiseProduct(b.At(x, y));
        }
    }
}

void Subtract(const Image3 &a, const Image3 &b, Image3 &out) {
    const int pixelWidth = a.pixelWidth;
    const int pixelHeight = a.pixelHeight;
    for (int y = 0; y < pixelHeight; y++) {
        for (int x = 0; x < pixelWidth; x++) {
            out.At(x, y) = a.At(x, y) - b.At(x, y);
        }
    }
}

void SubtractWithMultiplier(const Image3 &a, const Vector3 &multiplier, const Image3 &b, Image3 &out) {
    const int pixelWidth = a.pixelWidth;
    const int pixelHeight = a.pixelHeight;
    for (int y = 0; y < pixelHeight; y++) {
        for (int x = 0; x < pixelWidth; x++) {
            out.At(x, y) = a.At(x, y) - multiplier.cwiseProduct(b.At(x, y));
        }
    }
}

Vector3 xDotX(const Image3 &x_) {
    const int pixelWidth = x_.pixelWidth;
    const int pixelHeight = x_.pixelHeight;
    Vector3 result = Vector3::Zero();
    for (int y = 0; y < pixelHeight; y++) {
        for (int x = 0; x < pixelWidth; x++) {
            result += x_.At(x, y).cwiseProduct(x_.At(x, y));
        }
    }
    return result;
}

Vector3 xDotY(const Image3 &x_, const Image3 &y_) {
    const int pixelWidth = x_.pixelWidth;
    const int pixelHeight = x_.pixelHeight;
    Vector3 result = Vector3::Zero();
    for (int y = 0; y < pixelHeight; y++) {
        for (int x = 0; x < pixelWidth; x++) {
            result += x_.At(x, y).cwiseProduct(y_.At(x, y));
        }
    }
    return result;
}

void ComputeIRLSWeight(const double alpha,
                        const double irlsReg,
                        const Image3 &current,
                        const Image3 &throughputImage,
                        const Image3 &dxImage,
                        const Image3 &dyImage,
                        Image3 &weight) {
    const int pixelWidth = throughputImage.pixelWidth;
    const int pixelHeight = throughputImage.pixelHeight;
    for (int y = 0; y < pixelHeight; y++) {
        for (int x = 0; x < pixelWidth; x++) {
            const double throughputResidual =
                std::sqrt(alpha) * (throughputImage.At(x, y) - current.At(x, y)).norm();
            const Vector3 currentDx = x < pixelWidth - 1 ? 
                Vector3(current.At(x + 1, y) - current.At(x, y)) : 
                Vector3(-current.At(x, y));
            const double dxResidual =
                (dxImage.At(x, y) - currentDx).norm();
            const Vector3 currentDy = y < pixelHeight - 1 ? 
                Vector3(current.At(x, y + 1) - current.At(x, y)) :
                Vector3(- current.At(x, y));
            const double dyResidual =
                (dyImage.At(x, y) - currentDy).norm();
            weight.At(x, y)[0] = 1.0 / (throughputResidual + irlsReg);
            weight.At(x, y)[1] = 1.0 / (dxResidual + irlsReg);
            weight.At(x, y)[2] = 1.0 / (dyResidual + irlsReg);

        }
    }
}

void Solve(const Image3 &throughputImage,
           const Image3 &dxImage,
           const Image3 &dyImage,
           const Image3 &visibleEmitterImage,
           const double alpha,
           const bool doL1,
           Image3 &reconstructImage) {
    const int pixelWidth = throughputImage.pixelWidth;
    const int pixelHeight = throughputImage.pixelHeight;

    const double irlsRegInit = 0.05;
    const double irlsRegIter = 0.5;
    const int maxIrlsIteration = doL1 ? 20 : 1;

    Image3 x = throughputImage;
    Image3 throughputWeight(pixelWidth, pixelHeight);
    Image3 dxWeight(pixelWidth, pixelHeight);
    Image3 dyWeight(pixelWidth, pixelHeight);
    Image3 irlsWeight(pixelWidth, pixelHeight, Vector3::Constant(1.0));
    for (int irlsIteration = 0; irlsIteration < maxIrlsIteration; irlsIteration++) {
        if (irlsIteration > 0) {
            const double irlsReg = irlsRegInit * std::pow(irlsRegIter, (double)(irlsIteration - 1));
            ComputeIRLSWeight(alpha, irlsReg, x, throughputImage, dxImage, dyImage, irlsWeight);
        }

        Image3 atb(pixelWidth, pixelHeight);
        ComputeATb(alpha, throughputImage,
                   dxImage, dyImage, irlsWeight, atb);

        Image3 ataX = x;
        ApplyATA(alpha, x, irlsWeight, ataX);
        Image3 r = atb;
        Subtract(atb, ataX, r);
        Image3 p = r;
        Image3 ataP = p;
        for (int cgIteration = 0; cgIteration < 50; cgIteration++) {
            Vector3 rDotR = xDotX(r);
            if (rDotR.maxCoeff() <= 1e-20) {
                break;
            }
            ApplyATA(alpha, p, irlsWeight, ataP);
            Vector3 pDotAtaP = xDotY(p, ataP);
            Vector3 learningRate = rDotR.cwiseQuotient(pDotAtaP);
            AddWithMultiplier(x, learningRate, p, x);
            SubtractWithMultiplier(r, learningRate, ataP, r);
            Vector3 newRDotR = xDotX(r);
            Vector3 correctionFactor = newRDotR.cwiseQuotient(rDotR);
            AddWithMultiplier(r, correctionFactor, p, p);
        }
    }

    reconstructImage = x;

    for (int y = 0; y < pixelHeight; y++) {
        for (int x = 0; x < pixelWidth; x++) {
            reconstructImage.At(x, y) += visibleEmitterImage.At(x, y);
        }
    }
}

}
