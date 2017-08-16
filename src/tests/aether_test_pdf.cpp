#include <mitsuba/render/testcase.h>
#include <aether/RandomVar.h>
#include <aether/RandomSequence.h>
#include <aether/UniformTriangle.h>

using namespace aether;

std::atomic<int> aether::Object::next_id{0};

MTS_NAMESPACE_BEGIN

const Real c_testThreshold = 5e-3;

std::mt19937 rng{(uint32_t)time(NULL)};

template<typename Wi, typename M>
auto reflect(const Wi &wi, const M &m) {
    return two * dot(wi, m) * m - wi;
}

template<typename Wi, typename M>
auto refract(const Real eta, const Wi &wi, const M &m) {
    // http://www.graphics.cornell.edu/~bjw/microfacetbsdf.pdf
    // Eq 40, note that there is a typo and the eta inside the square root should be squared
    auto c = dot(wi, m);

    // eta = eta_interior / eta_exterior
    // eta_ = eta_incoming / eta_transmittance
    // if wi[2] > 0, we are at exterior, so eta_ = 1/eta, otherwise eta_ = eta
    auto eta_ = constant(at<2>(wi).Value() > Real(0) ? (Real(1) / eta) : eta);

    auto dsq = make_random_var(one) + sq(eta_) * (sq(c) - make_random_var(one));
    auto d = sqrt(dsq);
    auto e = pattern(when(at<2>(wi) > make_random_var(zero), d), otherwise(-d));
    auto ret = (eta_ * c - e) * m - eta_ * wi;
    return ret;
}

struct MyUniDist {
    Real Uniform1D() {
        return uniform(rng);
    }

    std::array<Real, 2> Uniform2D() {
        return std::array<Real, 2>{{uniform(rng), uniform(rng)}};
    }

    std::uniform_real_distribution<float> uniform{0.f, 1.f};
};

bool Compare(const std::string &name, const Real output, const Real expected) {
    bool is_same = fabs(output - expected) / (fabs(expected) + c_testThreshold) < c_testThreshold;
    if (!is_same) {
        std::cerr << "[Error] function:" << name << ", ouput:" << output << ", expected:" << expected << std::endl;
    }
    return is_same;
}

class TestAetherPdf : public TestCase {
public:
	MTS_BEGIN_TESTCASE()
	MTS_DECLARE_TEST(test_hemisphere)
	MTS_DECLARE_TEST(test_disk)
	MTS_DECLARE_TEST(test_triangle)
	MTS_DECLARE_TEST(test_spherical)
	MTS_DECLARE_TEST(test_cosine_hemisphere)
	MTS_DECLARE_TEST(test_phong)
    MTS_DECLARE_TEST(test_blinn_phong)
    MTS_DECLARE_TEST(test_cone)
    MTS_DECLARE_TEST(test_beckmann)
    MTS_DECLARE_TEST(test_beckmann_reflect)
    MTS_DECLARE_TEST(test_beckmann_refract)
    MTS_DECLARE_TEST(test_branching_brdf)
	MTS_END_TESTCASE()

	void test_hemisphere();
	void test_disk();
	void test_triangle();
	void test_spherical();
	void test_cosine_hemisphere();
	void test_phong();
    void test_blinn_phong();
    void test_cone();
    void test_beckmann();
    void test_beckmann_reflect();
    void test_beckmann_refract();
    void test_branching_brdf();
};

struct hemisphere_sampling_t {
    constexpr auto RV() const {
        constexpr auto r = sqrt(one - sq(u1));
        constexpr auto phi = two * get_pi() * u2;
        constexpr auto x = r * cos(phi);
        constexpr auto y = r * sin(phi);
        constexpr auto z = u1;
        return make_random_vector(make_random_var(x), make_random_var(y), make_random_var(z));
    }

    template <typename T>
    auto operator()(Context<T>& context, const RandomSequence<Vertex>& seq, MyUniDist& uniDist) const {
        auto uv = context.Uniform2D(uniDist);
        auto dir = RV().Sample(uv[0], uv[1]);

        return sample(dir);
    }
};

void TestAetherPdf::test_hemisphere() {
	for (int i = 0; i < 5; i++) {
		Node<hemisphere_sampling_t> hemisphere_sampling;
		RandomSequence<Vertex> path;
		MyUniDist uniDist;
		path.Append(hemisphere_sampling, uniDist);
		path.Sample();
		assertTrue(Compare(__func__, path.Pdf(), (1.f / (2.f * M_PI))));
	}
}

struct disk_sampling_t {
    constexpr auto RV() const {
        constexpr auto r = sqrt(u1);
        constexpr auto phi = two * get_pi() * u2;
        constexpr auto x = r * cos(phi);
        constexpr auto y = r * sin(phi);
        constexpr auto z = zero;
        return make_random_vector(make_random_var(x), make_random_var(y), make_random_var(z));
    }

    template <typename T>
    auto operator()(Context<T>& context, const RandomSequence<Vertex>& seq, MyUniDist& uniDist) const {
        auto uv = context.Uniform2D(uniDist);
        auto pos = RV().Sample(uv[0], uv[1]);

        return sample(pos);
    }
};

void TestAetherPdf::test_disk() {
	for (int i = 0; i < 5; i++) {
		Node<disk_sampling_t> disk_sampling;
		RandomSequence<Vertex> path;
		MyUniDist uniDist;
		path.Append(disk_sampling, uniDist);
		path.Sample();
		assertTrue(Compare(__func__, path.Pdf(), (1.f / M_PI)));
	}
}

struct triangle_sampling_t {
    template <typename T>
    auto operator()(Context<T>& context, const RandomSequence<Vertex>& seq, MyUniDist& uniDist) const {
        auto uv = context.Uniform2D(uniDist);
        auto v0 = aether::Vector3(0.0, 0.0, 0.0);
        auto v1 = aether::Vector3(1.0, 0.0, 0.0);
        auto v2 = aether::Vector3(0.0, 1.0, 0.0);
        auto pos = uniform_triangle(v0, v1, v2).Sample(uv[0], uv[1]);

        return sample(pos);
    }
};

void TestAetherPdf::test_triangle() {
	for (int i = 0; i < 5; i++) {
		Node<triangle_sampling_t> triangle_sampling;
		RandomSequence<Vertex> path;
		MyUniDist uniDist;
		path.Append(triangle_sampling, uniDist);
		path.Sample();
		assertTrue(Compare(__func__, path.Pdf(), 1.f / 0.5f));
	}
}

struct spherical_sampling_t {
    constexpr auto RV() const {
        constexpr auto phi = two * get_pi() * u1;
        constexpr auto sinPhi = sin(phi);
        constexpr auto cosPhi = cos(phi);
        constexpr auto theta = get_pi() * u2;
        constexpr auto sinTheta = sin(theta);
        constexpr auto cosTheta = cos(theta);
        return make_random_vector(
            make_random_var(sinPhi * sinTheta),
            make_random_var(cosTheta),
            make_random_var(- cosPhi * sinTheta));
    }

    template <typename T>
    auto operator()(Context<T>& context, const RandomSequence<Vertex>& seq, MyUniDist& uniDist) const {
        auto uv = context.Uniform2D(uniDist);
        auto dir = RV().Sample(uv[0], uv[1]);

        return sample(dir);
    }
};

void TestAetherPdf::test_spherical() {
	for (int i = 0; i < 5; i++) {
        Node<spherical_sampling_t> spherical_sampling;
		RandomSequence<Vertex> path;
		MyUniDist uniDist;
		path.Append(spherical_sampling, uniDist);
		path.Sample();
		// We sample uniformly in spherical coordinates, pdf is 1.0 / (2pi * pi * sin(theta))
		const Real theta = acos(path.Back().Value()[1]);
		assertTrue(Compare(__func__, path.Pdf(), fabs(1.0 / (2.f * M_PI * M_PI * sin(theta)))));
	}
}

struct cosine_hemisphere_sampling_t {
    constexpr auto RV() const {
        constexpr auto r = sqrt(one - u1);
        constexpr auto phi = two * get_pi() * u2;
        constexpr auto x = r * cos(phi);
        constexpr auto y = r * sin(phi);
        constexpr auto z = sqrt(u1);
        return make_random_vector(make_random_var(x), make_random_var(y), make_random_var(z));
    }

    template <typename T>
    auto operator()(Context<T>& context, const RandomSequence<Vertex>& seq, MyUniDist& uniDist) const {
        auto uv = context.Uniform2D(uniDist);
        auto dir = RV().Sample(uv[0], uv[1]);

        return sample(dir);
    }
};

void TestAetherPdf::test_cosine_hemisphere() {
	for (int i = 0; i < 5; i++) {
		Node<cosine_hemisphere_sampling_t> cosine_hemisphere_sampling;
		RandomSequence<Vertex> path;
		MyUniDist uniDist;
		path.Append(cosine_hemisphere_sampling, uniDist);
		path.Sample();
		assertTrue(Compare(__func__, path.Pdf(), fabs(path.Back().Value()[2] / M_PI)));
	}
}

struct phong_sampling_t {
    auto RV() const {
        auto cosine = pow(make_random_var(u1), constant(1.0 / (power + 1.0)));
        auto r = sqrt(one - sq(cosine));
        auto phi = two * get_pi() * u2;
        auto x = r * cos(phi);
        auto y = r * sin(phi);
        auto z = cosine;
        return make_random_vector(x, y, z);
    }

    template <typename T>
    auto operator()(Context<T>& context, const RandomSequence<Vertex>& seq, MyUniDist& uniDist) const {
        auto uv = context.Uniform2D(uniDist);
        auto dir = RV().Sample(uv[0], uv[1]);

        return sample(dir);
    }

    Real power;
};

void TestAetherPdf::test_phong() {
    aether::Real power[] = {1.0, 2.0, 4.0, 8.0, 12.0};
    for (int i = 0; i < 5; i++) {
		Node<phong_sampling_t> phong_sampling{power[i]};
		RandomSequence<Vertex> path;
		MyUniDist uniDist;
		path.Append(phong_sampling, uniDist);
		path.Sample();
		const Real z = path.Back().Value()[2];
		const Real analyticalPdf = (power[i] + 1.f) * powf(z, power[i]) * (0.5f / M_PI);
		assertTrue(Compare(__func__, path.Pdf(), analyticalPdf));
	}
}

struct blinn_phong_sampling_t {
    auto RV() const {
        auto cosine = pow(make_random_var(u1), constant(1.0 / (power + 1.0)));
        auto r = sqrt(one - sq(cosine));
        auto phi = two * get_pi() * u2;
        auto x = r * cos(phi);
        auto y = r * sin(phi);
        auto z = cosine;
        return make_random_vector(x, y, z);
    }

    template <typename T>
    auto operator()(Context<T>& context, const RandomSequence<Vertex>& seq, MyUniDist& uniDist) const {
        auto uv = context.Uniform2D(uniDist);
        auto wi = constant(wi_v);
        auto dir = reflect(wi, RV().Sample(uv[0], uv[1]));

        return sample(dir);
    }

    aether::Vector3 wi_v;
    Real power;
};

void TestAetherPdf::test_blinn_phong() {
    aether::Real power[] = {1.0, 2.0, 4.0, 8.0, 12.0};
    aether::Vector3 wi = aether::Vector3(0.1, 0.2, 0.9).normalized();
    for (int i = 0; i < 5; i++) {
        Node<blinn_phong_sampling_t> blinn_phong_sampling{wi, power[i]};
        RandomSequence<Vertex> path;
        MyUniDist uniDist;
        path.Append(blinn_phong_sampling, uniDist);
        path.Sample();
        const aether::Vector3 H = (wi + path.Back().Value()).normalized();
        const Real dwh_dwo = 1.0f / (4.0f * path.Back().Value().dot(H));
        const Real cosTheta = fabs(H[2]);
        const Real analyticalPdf = dwh_dwo * (power[i] + 1.f) * powf(cosTheta, power[i]) * (0.5f / M_PI);
        assertTrue(Compare(__func__, path.Pdf(), analyticalPdf));
    }
}

struct cone_sampling_t {
    auto RV() const {
        auto u1r = make_random_var(u1);
        auto cosTheta = make_random_var(one) + u1r * constant(cosCutOff - 1.0);
        auto sinTheta = sqrt(one - sq(cosTheta));
        auto phi = two * get_pi() * u2;
        auto sinPhi = sin(phi);
        auto cosPhi = cos(phi);
        return make_random_vector(
            cosPhi * sinTheta,
            sinPhi * sinTheta,
            cosTheta
        );
    }

    template <typename T>
    auto operator()(Context<T>& context, const RandomSequence<Vertex>& seq, MyUniDist& uniDist) const {
        auto uv = context.Uniform2D(uniDist);
        auto dir = RV().Sample(uv[0], uv[1]);

        return sample(dir);
    }

    Real cosCutOff;
};

void TestAetherPdf::test_cone() {
    aether::Real cosCutOff[] = {0.1, 0.2, 0.4, 0.6, 0.8};
    for (int i = 0; i < 5; i++) {
        Node<cone_sampling_t> cone_sampling{cosCutOff[i]};
        RandomSequence<Vertex> path;
        MyUniDist uniDist;
        path.Append(cone_sampling, uniDist);
        path.Sample();
        const Real analyticalPdf = 1.0 / (2.0 * M_PI * (1.0 - cosCutOff[i]));
        assertTrue(Compare(__func__, path.Pdf(), analyticalPdf));
    }
}

struct beckmann_sampling_t {
    auto RV(const aether::Real alpha) const {
        auto phi = two * get_pi() * u2;
        auto tanThetaSqr = constant(- alpha * alpha) * make_random_var(log(one - u1));
        auto cosTheta = make_random_var(one) / sqrt(make_random_var(one) + tanThetaSqr);
        auto sinTheta = sqrt(make_random_var(one) - sq(cosTheta));
        return make_random_vector(
            sinTheta * cos(phi),
            sinTheta * sin(phi),
            cosTheta
        );
    }

    template <typename T>
    auto operator()(Context<T>& context, const RandomSequence<Vertex>& seq, MyUniDist& uniDist) const {
        auto uv = context.Uniform2D(uniDist);
        auto dir = RV(alpha).Sample(uv[0], uv[1]);

        return sample(dir);
    }

    Real alpha;
};

void TestAetherPdf::test_beckmann() {
    aether::Real alpha[] = {0.01, 0.05, 0.1, 0.2, 0.4};
    for (int i = 0; i < 5; i++) {
        Node<beckmann_sampling_t> beckmann_sampling{alpha[i]};
        RandomSequence<Vertex> path;
        MyUniDist uniDist;
        path.Append(beckmann_sampling, uniDist);
        path.Sample();
        const Real cosTheta = path.Back().Value()[2];
        const Real theta = acos(cosTheta);
        const Real tanTheta = tan(theta);
        const Real exponent = exp(-(tanTheta * tanTheta) / (alpha[i] * alpha[i]));
        const Real analyticalPdf = exponent / (M_PI * alpha[i] * alpha[i] * cosTheta * cosTheta * cosTheta);
        assertTrue(Compare(__func__, path.Pdf(), analyticalPdf));
    }
}

struct beckmann_sampling_reflect_t {
    auto RV(const aether::Real alpha) const {
        auto phi = two * get_pi() * u2;
        auto tanThetaSqr = constant(- alpha * alpha) * make_random_var(log(one - u1));
        auto cosTheta = make_random_var(one) / sqrt(make_random_var(one) + tanThetaSqr);
        auto sinTheta = sqrt(make_random_var(one) - sq(cosTheta));
        return make_random_vector(
            sinTheta * cos(phi),
            sinTheta * sin(phi),
            cosTheta
        );
    }

    template <typename T>
    auto operator()(Context<T>& context, const RandomSequence<Vertex>& seq, MyUniDist& uniDist) const {
        auto uv = context.Uniform2D(uniDist);
        auto wi = constant(wi_v);
        auto dir = reflect(wi, RV(alpha).Sample(uv[0], uv[1]));

        return sample(dir);
    }

    aether::Vector3 wi_v;
    Real alpha;
};

void TestAetherPdf::test_beckmann_reflect() {
    aether::Real alpha[] = {0.01, 0.05, 0.1, 0.2, 0.4};
    aether::Vector3 wi = aether::Vector3(0.1, 0.2, 0.9).normalized();
    for (int i = 0; i < 5; i++) {
        Node<beckmann_sampling_reflect_t> beckmann_sampling{wi, alpha[i]};
        RandomSequence<Vertex> path;
        MyUniDist uniDist;
        path.Append(beckmann_sampling, uniDist);
        path.Sample();
        const aether::Vector3 H = (wi + path.Back().Value()).normalized();
        const Real dwh_dwo = 1.0f / (4.0f * path.Back().Value().dot(H));
        const Real cosTheta = H[2];
        const Real theta = acos(cosTheta);
        const Real tanTheta = tan(theta);
        const Real exponent = exp(-(tanTheta * tanTheta) / (alpha[i] * alpha[i]));
        const Real analyticalPdf = dwh_dwo * exponent / (M_PI * alpha[i] * alpha[i] * cosTheta * cosTheta * cosTheta);
        assertTrue(Compare(__func__, path.Pdf(), analyticalPdf));
    }
}

struct beckmann_sampling_refract_t {
    auto RV(const aether::Real alpha) const {
        auto phi = two * get_pi() * u2;
        auto tanThetaSqr = constant(- alpha * alpha) * make_random_var(log(one - u1));
        auto cosTheta = make_random_var(one) / sqrt(make_random_var(one) + tanThetaSqr);
        auto sinTheta = sqrt(make_random_var(one) - sq(cosTheta));
        return make_random_vector(
            sinTheta * cos(phi),
            sinTheta * sin(phi),
            cosTheta
        );
    }

    template <typename T>
    auto operator()(Context<T>& context, const RandomSequence<Vertex>& seq, MyUniDist& uniDist) const {
        auto uv = context.Uniform2D(uniDist);
        auto wi = constant(wi_v);
        auto dir = refract(eta, wi, RV(alpha).Sample(uv[0], uv[1]));

        return sample(dir);
    }

    Real eta;
    aether::Vector3 wi_v;
    Real alpha;
};

void TestAetherPdf::test_beckmann_refract() {
    aether::Real alpha[] = {0.01, 0.05, 0.1, 0.2, 0.4};
    aether::Vector3 wi = aether::Vector3(0.1, 0.2, 0.9).normalized();
    for (int i = 0; i < 5; i++) {
        std::uniform_real_distribution<float> uniform{0.f, 1.f};
        aether::Real eta = uniform(rng) + 0.5f;
        Node<beckmann_sampling_refract_t> beckmann_sampling{eta, wi, alpha[i]};
        RandomSequence<Vertex> path;
        MyUniDist uniDist;
        path.Append(beckmann_sampling, uniDist);
        path.Sample();
        const aether::Vector3 H = (wi + path.Back().Value() * eta).normalized();
        const Real sqrtDenom = wi.dot(H) + eta * path.Back().Value().dot(H);
        const Real dwh_dwo = fabs((eta*eta * path.Back().Value().dot(H)) / (sqrtDenom*sqrtDenom));
        const Real cosTheta = fabs(H[2]);
        const Real theta = acos(cosTheta);
        const Real tanTheta = tan(theta);
        const Real exponent = exp(-(tanTheta * tanTheta) / (alpha[i] * alpha[i]));
        const Real analyticalPdf = dwh_dwo * exponent / (M_PI * alpha[i] * alpha[i] * cosTheta * cosTheta * cosTheta);
        assertTrue(Compare(__func__, path.Pdf(), analyticalPdf));
    }
}

template<typename S1, typename S2>
struct branching_sampling_t {
    template <typename T>
    auto operator()(Context<T>& context, const RandomSequence<Vertex>& seq, MyUniDist& uniDist) const {
        auto componentSampler = discrete_dynamic(std::vector<Real>{0, 1});
        auto component = context.Sample(componentSampler, context.Uniform1D(uniDist));
        auto dir = pattern(
            when(constant(component) == make_random_var(zero), sampler1(context, seq, uniDist)),
            otherwise(sampler2(context, seq, uniDist))
        );

        return sample(dir);
    }

    S1 sampler1;
    S2 sampler2;
};

void TestAetherPdf::test_branching_brdf() {
    // uniform hemisphere + uniform hemisphere
    for (int i = 0; i < 5; i++) {
        Node<branching_sampling_t<hemisphere_sampling_t, hemisphere_sampling_t>> 
            branching_sampling{hemisphere_sampling_t{}, hemisphere_sampling_t{}};
        RandomSequence<Vertex> path;
        MyUniDist uniDist;
        path.Append(branching_sampling, uniDist);
        path.Sample();
        const Real hemispherePdf = fabs(1.0 / (2.f * M_PI));
        assertTrue(Compare(__func__, path.Pdf(), hemispherePdf));
    }

    // cosine hemisphere + uniform hemisphere
    for (int i = 0; i < 5; i++) {
        Node<branching_sampling_t<cosine_hemisphere_sampling_t, hemisphere_sampling_t>> 
            branching_sampling{cosine_hemisphere_sampling_t{}, hemisphere_sampling_t{}};
        RandomSequence<Vertex> path;
        MyUniDist uniDist;
        path.Append(branching_sampling, uniDist);
        path.Sample();
        const Real hemispherePdf = fabs(1.0 / (2.f * M_PI));
        const Real coshemispherePdf = fabs(path.Back().Value()[2] / M_PI);
        assertTrue(Compare(__func__, path.Pdf(), (hemispherePdf + coshemispherePdf) * Real(0.5)));
    }

    aether::Real alpha[] = {0.01, 0.05, 0.1, 0.2, 0.4};
    aether::Real power[] = {1.0, 2.0, 4.0, 8.0, 12.0};
    aether::Vector3 wi = aether::Vector3(0.1, 0.2, 0.9).normalized();
    // uniform hemisphere + phong sampling
    for (int i = 0; i < 5; i++) {
        Node<branching_sampling_t<hemisphere_sampling_t, phong_sampling_t>> 
            branching_sampling{hemisphere_sampling_t{}, phong_sampling_t{power[i]}};
        RandomSequence<Vertex> path;
        MyUniDist uniDist;
        path.Append(branching_sampling, uniDist);
        path.Sample();
        const Real z = path.Back().Value()[2];
        const Real phongPdf = (power[i] + 1.f) * powf(z, power[i]) * (0.5f / M_PI);
        const Real hemispherePdf = fabs(1.0 / (2.f * M_PI));
        assertTrue(Compare(__func__, path.Pdf(), (phongPdf + hemispherePdf) * Real(0.5)));
    }

    // cosine uniform hemisphere + blinn-phong sampling
    for (int i = 0; i < 5; i++) {
        Node<branching_sampling_t<cosine_hemisphere_sampling_t, blinn_phong_sampling_t>> 
            branching_sampling{cosine_hemisphere_sampling_t{}, blinn_phong_sampling_t{wi, power[i]}};
        RandomSequence<Vertex> path;
        MyUniDist uniDist;
        path.Append(branching_sampling, uniDist);
        path.Sample();
        const aether::Vector3 H = (wi + path.Back().Value()).normalized();
        const Real dwh_dwo = 1.0f / (4.0f * path.Back().Value().dot(H));
        const Real cosTheta = fabs(H[2]);
        const Real blinnPhongPdf = dwh_dwo * (power[i] + 1.f) * powf(cosTheta, power[i]) * (0.5f / M_PI);
        const Real coshemispherePdf = fabs(path.Back().Value()[2] / M_PI);
        assertTrue(Compare(__func__, path.Pdf(), (blinnPhongPdf + coshemispherePdf) * Real(0.5)));
    }
    
    // uniform hemisphere + beckmann sampling
    for (int i = 0; i < 5; i++) {
        Node<branching_sampling_t<hemisphere_sampling_t, beckmann_sampling_t>> 
            branching_sampling{hemisphere_sampling_t{}, beckmann_sampling_t{alpha[i]}};
        RandomSequence<Vertex> path;
        MyUniDist uniDist;
        path.Append(branching_sampling, uniDist);
        path.Sample();
        const Real cosTheta = path.Back().Value()[2];
        const Real theta = acos(cosTheta);
        const Real tanTheta = tan(theta);
        const Real exponent = exp(-(tanTheta * tanTheta) / (alpha[i] * alpha[i]));
        const Real beckmannPdf = exponent / (M_PI * alpha[i] * alpha[i] * cosTheta * cosTheta * cosTheta);
        const Real hemispherePdf = fabs(1.0 / (2.f * M_PI));
        assertTrue(Compare(__func__, path.Pdf(), (beckmannPdf + hemispherePdf) * Real(0.5)));
    }

    // cosine hemisphere + beckmann sampling
    for (int i = 0; i < 5; i++) {
        Node<branching_sampling_t<cosine_hemisphere_sampling_t, beckmann_sampling_reflect_t>> 
            branching_sampling{cosine_hemisphere_sampling_t{}, beckmann_sampling_reflect_t{wi, alpha[i]}};
        RandomSequence<Vertex> path;
        MyUniDist uniDist;
        path.Append(branching_sampling, uniDist);
        path.Sample();
        const aether::Vector3 H = (wi + path.Back().Value()).normalized();
        const Real dwh_dwo = 1.0f / (4.0f * path.Back().Value().dot(H));
        const Real cosTheta = H[2];
        const Real theta = acos(cosTheta);
        const Real tanTheta = tan(theta);
        const Real exponent = exp(-(tanTheta * tanTheta) / (alpha[i] * alpha[i]));
        const Real beckmannAnalyticalPdf = dwh_dwo * exponent / (M_PI * alpha[i] * alpha[i] * cosTheta * cosTheta * cosTheta);
        const Real coshemispherePdf = fabs(path.Back().Value()[2] / M_PI);
        assertTrue(Compare(__func__, path.Pdf(), (beckmannAnalyticalPdf + coshemispherePdf) * Real(0.5)));
    }
}

MTS_EXPORT_TESTCASE(TestAetherPdf, "Testcase for aether pdf")

MTS_NAMESPACE_END
