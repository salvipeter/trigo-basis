#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

#include <QtGui/QKeyEvent>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Tools/Smoother/JacobiLaplaceSmootherT.hh>

// #define BETTER_MEAN_CURVATURE

#ifdef BETTER_MEAN_CURVATURE
#include "Eigen/Eigenvalues"
#include "Eigen/Geometry"
#include "Eigen/LU"
#include "Eigen/SVD"
#endif

#include "MyViewer.h"

#ifdef _WIN32
#define GL_CLAMP_TO_EDGE 0x812F
#define GL_BGRA 0x80E1
#endif

MyViewer::MyViewer(QWidget *parent) :
  QGLViewer(parent), model_type(ModelType::NONE),
  trigonometric_basis(false),
  mean_min(0.0), mean_max(0.0), cutoff_ratio(0.05),
  show_control_points(true), show_solid(true), show_wireframe(false),
  visualization(Visualization::PLAIN), slicing_dir(0, 0, 1), slicing_scaling(1),
  last_filename("")
{
  setSelectRegionWidth(10);
  setSelectRegionHeight(10);
  axes.shown = false;
}

MyViewer::~MyViewer() {
  glDeleteTextures(1, &isophote_texture);
  glDeleteTextures(1, &environment_texture);
  glDeleteTextures(1, &slicing_texture);
}

void MyViewer::updateMeanMinMax() {
  size_t n = mesh.n_vertices();
  if (n == 0)
    return;

  std::vector<double> mean;
  mean.reserve(n);
  for (auto v : mesh.vertices())
    mean.push_back(mesh.data(v).mean);

  std::sort(mean.begin(), mean.end());
  size_t k = (double)n * cutoff_ratio;
  mean_min = std::min(mean[k ? k-1 : 0], 0.0);
  mean_max = std::max(mean[n-k], 0.0);
}

void MyViewer::localSystem(const MyViewer::Vector &normal,
                           MyViewer::Vector &u, MyViewer::Vector &v) {
  // Generates an orthogonal (u,v) coordinate system in the plane defined by `normal`.
  int maxi = 0, nexti = 1;
  double max = std::abs(normal[0]), next = std::abs(normal[1]);
  if (max < next) {
    std::swap(max, next);
    std::swap(maxi, nexti);
  }
  if (std::abs(normal[2]) > max) {
    nexti = maxi;
    maxi = 2;
  } else if (std::abs(normal[2]) > next)
    nexti = 2;

  u.vectorize(0.0);
  u[nexti] = -normal[maxi];
  u[maxi] = normal[nexti];
  u /= u.norm();
  v = normal % u;
}

double MyViewer::voronoiWeight(MyViewer::MyMesh::HalfedgeHandle in_he) {
  // Returns the area of the triangle bounded by in_he that is closest
  // to the vertex pointed to by in_he.
  if (mesh.is_boundary(in_he))
    return 0;
  auto next = mesh.next_halfedge_handle(in_he);
  auto prev = mesh.prev_halfedge_handle(in_he);
  double c2 = mesh.calc_edge_vector(in_he).sqrnorm();
  double b2 = mesh.calc_edge_vector(next).sqrnorm();
  double a2 = mesh.calc_edge_vector(prev).sqrnorm();
  double alpha = mesh.calc_sector_angle(in_he);

  if (a2 + b2 < c2)                // obtuse gamma
    return 0.125 * b2 * std::tan(alpha);
  if (a2 + c2 < b2)                // obtuse beta
    return 0.125 * c2 * std::tan(alpha);
  if (b2 + c2 < a2) {              // obtuse alpha
    double b = std::sqrt(b2), c = std::sqrt(c2);
    double total_area = 0.5 * b * c * std::sin(alpha);
    double beta  = mesh.calc_sector_angle(prev);
    double gamma = mesh.calc_sector_angle(next);
    return total_area - 0.125 * (b2 * std::tan(gamma) + c2 * std::tan(beta));
  }

  double r2 = 0.25 * a2 / std::pow(std::sin(alpha), 2); // squared circumradius
  auto area = [r2](double x2) {
    return 0.125 * std::sqrt(x2) * std::sqrt(std::max(4.0 * r2 - x2, 0.0));
  };
  return area(b2) + area(c2);
}

#ifndef BETTER_MEAN_CURVATURE
void MyViewer::updateMeanCurvature(bool update_min_max) {
  if (model_type == ModelType::BEZIER_SURFACE) {
    for (auto v : mesh.vertices()) {
      std::vector<std::vector<Vec>> der;
      evaluate(mesh.data(v).u, mesh.data(v).v, 2, der);
      auto E = der[1][0].squaredNorm();
      auto F = der[1][0] * der[0][1];
      auto G = der[0][1].squaredNorm();
      auto n = (der[1][0] ^ der[0][1]).unit();
      auto L = n * der[2][0];
      auto M = n * der[1][1];
      auto N = n * der[0][2];
      mesh.data(v).mean = (N * E - 2 * M * F + L * G) / (2 * (E * G - F * F));
      // mesh.data(v).gauss = (L * N - M * M) / (E * G - F * F);
    }

    if (update_min_max)
      updateMeanMinMax();
    return;
  }

  std::map<MyMesh::FaceHandle, double> face_area;
  std::map<MyMesh::VertexHandle, double> vertex_area;

  for (auto f : mesh.faces())
    face_area[f] = mesh.calc_sector_area(mesh.halfedge_handle(f));

  // Compute triangle strip areas
  for (auto v : mesh.vertices()) {
    vertex_area[v] = 0;
    mesh.data(v).mean = 0;
    for (auto f : mesh.vf_range(v))
      vertex_area[v] += face_area[f];
    vertex_area[v] /= 3.0;
  }

  // Compute mean values using dihedral angles
  for (auto v : mesh.vertices()) {
    for (auto h : mesh.vih_range(v)) {
      auto vec = mesh.calc_edge_vector(h);
      double angle = mesh.calc_dihedral_angle(h); // signed; returns 0 at the boundary
      mesh.data(v).mean += angle * vec.norm();
    }
    mesh.data(v).mean *= 0.25 / vertex_area[v];
  }

  if (update_min_max)
    updateMeanMinMax();
}
#else // BETTER_MEAN_CURVATURE
void MyViewer::updateMeanCurvature(bool update_min_max) {
  // As in the paper:
  //   S. Rusinkiewicz, Estimating curvatures and their derivatives on triangle meshes.
  //     3D Data Processing, Visualization and Transmission, IEEE, 2004.

  std::map<MyMesh::VertexHandle, Vector> efgp; // 2nd principal form
  std::map<MyMesh::VertexHandle, double> wp;   // accumulated weight

  // Initial setup
  for (auto v : mesh.vertices()) {
    efgp[v].vectorize(0.0);
    wp[v] = 0.0;
  }

  for (auto f : mesh.faces()) {
    // Setup local edges, vertices and normals
    auto h0 = mesh.halfedge_handle(f);
    auto h1 = mesh.next_halfedge_handle(h0);
    auto h2 = mesh.next_halfedge_handle(h1);
    auto e0 = mesh.calc_edge_vector(h0);
    auto e1 = mesh.calc_edge_vector(h1);
    auto e2 = mesh.calc_edge_vector(h2);
    auto n0 = mesh.normal(mesh.to_vertex_handle(h1));
    auto n1 = mesh.normal(mesh.to_vertex_handle(h2));
    auto n2 = mesh.normal(mesh.to_vertex_handle(h0));

    Vector n = mesh.normal(f), u, v;
    localSystem(n, u, v);

    // Solve a LSQ equation for (e,f,g) of the face
    Eigen::MatrixXd A(6, 3);
    A << (e0 | u), (e0 | v),    0.0,
            0.0,   (e0 | u), (e0 | v),
         (e1 | u), (e1 | v),    0.0,
            0.0,   (e1 | u), (e1 | v),
         (e2 | u), (e2 | v),    0.0,
            0.0,   (e2 | u), (e2 | v);
    Eigen::VectorXd b(6);
    b << ((n2 - n1) | u),
         ((n2 - n1) | v),
         ((n0 - n2) | u),
         ((n0 - n2) | v),
         ((n1 - n0) | u),
         ((n1 - n0) | v);
    Eigen::Vector3d x = A.fullPivLu().solve(b);

    Eigen::Matrix2d F;          // Fundamental matrix for the face
    F << x(0), x(1),
         x(1), x(2);

    for (auto h : mesh.fh_range(f)) {
      auto p = mesh.to_vertex_handle(h);

      // Rotate the (up,vp) local coordinate system to be coplanar with that of the face
      Vector np = mesh.normal(p), up, vp;
      localSystem(np, up, vp);
      auto axis = (np % n).normalize();
      double angle = std::acos(std::min(std::max(n | np, -1.0), 1.0));
      auto rotation = Eigen::AngleAxisd(angle, Eigen::Vector3d(axis.data()));
      Eigen::Vector3d up1(up.data()), vp1(vp.data());
      up1 = rotation * up1;    vp1 = rotation * vp1;
      up = Vector(up1.data()); vp = Vector(vp1.data());

      // Compute the vertex-local (e,f,g)
      double e, f, g;
      Eigen::Vector2d upf, vpf;
      upf << (up | u), (up | v);
      vpf << (vp | u), (vp | v);
      e = upf.transpose() * F * upf;
      f = upf.transpose() * F * vpf;
      g = vpf.transpose() * F * vpf;

      // Accumulate the results with Voronoi weights
      double w = voronoiWeight(h);
      efgp[p] += Vector(e, f, g) * w;
      wp[p] += w;
    }
  }

  // Compute the principal curvatures
  for (auto v : mesh.vertices()) {
    auto &efg = efgp[v];
    efg /= wp[v];
    Eigen::Matrix2d F;
    F << efg[0], efg[1],
         efg[1], efg[2];
    auto k = F.eigenvalues();   // always real, because F is a symmetric real matrix
    mesh.data(v).mean = (k(0).real() + k(1).real()) / 2.0;
  }

  if (update_min_max)
    updateMeanMinMax();
}
#endif

static Vec HSV2RGB(Vec hsv) {
  // As in Wikipedia
  double c = hsv[2] * hsv[1];
  double h = hsv[0] / 60;
  double x = c * (1 - std::abs(std::fmod(h, 2) - 1));
  double m = hsv[2] - c;
  Vec rgb(m, m, m);
  if (h <= 1)
    return rgb + Vec(c, x, 0);
  if (h <= 2)
    return rgb + Vec(x, c, 0);
  if (h <= 3)
    return rgb + Vec(0, c, x);
  if (h <= 4)
    return rgb + Vec(0, x, c);
  if (h <= 5)
    return rgb + Vec(x, 0, c);
  if (h <= 6)
    return rgb + Vec(c, 0, x);
  return rgb;
}

Vec MyViewer::meanMapColor(double d) const {
  double red = 0, green = 120, blue = 240; // Hue
  if (d < 0) {
    double alpha = mean_min ? std::min(d / mean_min, 1.0) : 1.0;
    return HSV2RGB({green * (1 - alpha) + blue * alpha, 1, 1});
  }
  double alpha = mean_max ? std::min(d / mean_max, 1.0) : 1.0;
  return HSV2RGB({green * (1 - alpha) + red * alpha, 1, 1});
}

void MyViewer::fairMesh() {
  if (model_type != ModelType::MESH)
    return;

  emit startComputation(tr("Fairing mesh..."));
  OpenMesh::Smoother::JacobiLaplaceSmootherT<MyMesh> smoother(mesh);
  smoother.initialize(OpenMesh::Smoother::SmootherT<MyMesh>::Normal, // or: Tangential_and_Normal
                      OpenMesh::Smoother::SmootherT<MyMesh>::C1);
  for (size_t i = 1; i <= 10; ++i) {
    smoother.smooth(10);
    emit midComputation(i * 10);
  }
  updateMesh(false);
  emit endComputation();
}

void MyViewer::updateVertexNormals() {
  if (model_type == ModelType::BEZIER_SURFACE) {
    for (auto v : mesh.vertices()) {
      std::vector<std::vector<Vec>> der;
      evaluate(mesh.data(v).u, mesh.data(v).v, 1, der);
      Vector n(static_cast<double *>(der[0][1] ^ der[1][0]));
      double len = n.length();
      if (len != 0.0)
        n /= len;
      mesh.set_normal(v, n);
    }
    return;
  }

  // Weights according to:
  //   N. Max, Weights for computing vertex normals from facet normals.
  //     Journal of Graphics Tools, Vol. 4(2), 1999.
  for (auto v : mesh.vertices()) {
    Vector n(0.0, 0.0, 0.0);
    for (auto h : mesh.vih_range(v)) {
      if (mesh.is_boundary(h))
        continue;
      auto in_vec  = mesh.calc_edge_vector(h);
      auto out_vec = mesh.calc_edge_vector(mesh.next_halfedge_handle(h));
      double w = in_vec.sqrnorm() * out_vec.sqrnorm();
      n += (in_vec % out_vec) / (w == 0.0 ? 1.0 : w);
    }
    double len = n.length();
    if (len != 0.0)
      n /= len;
    mesh.set_normal(v, n);
  }
}

void MyViewer::updateMesh(bool update_mean_range) {
  if (model_type == ModelType::BEZIER_SURFACE)
    generateMesh(50);
  mesh.request_face_normals(); mesh.request_vertex_normals();
  mesh.update_face_normals(); //mesh.update_vertex_normals();
  updateVertexNormals();
  updateMeanCurvature(update_mean_range);
}

void MyViewer::setupCamera() {
  // Set camera on the model
  Vector box_min, box_max;
  box_min = box_max = mesh.point(*mesh.vertices_begin());
  for (auto v : mesh.vertices()) {
    box_min.minimize(mesh.point(v));
    box_max.maximize(mesh.point(v));
  }
  camera()->setSceneBoundingBox(Vec(box_min.data()), Vec(box_max.data()));
  camera()->showEntireScene();

  slicing_scaling = 20 / (box_max - box_min).max();

  setSelectedName(-1);
  axes.shown = false;

  update();
}

bool MyViewer::openMesh(const std::string &filename, bool update_view) {
  if (!OpenMesh::IO::read_mesh(mesh, filename) || mesh.n_vertices() == 0)
    return false;
  model_type = ModelType::MESH;
  last_filename = filename;
  updateMesh(update_view);
  if (update_view)
    setupCamera();
  return true;
}

bool MyViewer::openBezier(const std::string &filename, bool update_view) {
  size_t n, m;
  try {
    std::ifstream f(filename.c_str());
    f.exceptions(std::ios::failbit | std::ios::badbit);
    f >> n >> m;
    degree[0] = n++; degree[1] = m++;
    control_points.resize(n * m);
    for (size_t i = 0, index = 0; i < n; ++i)
      for (size_t j = 0; j < m; ++j, ++index)
        f >> control_points[index][0] >> control_points[index][1] >> control_points[index][2];
  } catch(std::ifstream::failure &) {
    return false;
  }
  model_type = ModelType::BEZIER_SURFACE;
  last_filename = filename;
  updateMesh(update_view);
  if (update_view)
    setupCamera();
  return true;
}

bool MyViewer::saveBezier(const std::string &filename) {
  if (model_type != ModelType::BEZIER_SURFACE)
    return false;

  try {
    std::ofstream f(filename.c_str());
    f.exceptions(std::ios::failbit | std::ios::badbit);
    f << degree[0] << ' ' << degree[1] << std::endl;
    for (const auto &p : control_points)
      f << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
  } catch(std::ifstream::failure &) {
    return false;
  }
  return true;
}

void MyViewer::init() {
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);

  QImage img(":/isophotes.png");
  glGenTextures(1, &isophote_texture);
  glBindTexture(GL_TEXTURE_2D, isophote_texture);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, img.width(), img.height(), 0, GL_BGRA,
               GL_UNSIGNED_BYTE, img.convertToFormat(QImage::Format_ARGB32).bits());

  QImage img2(":/environment.png");
  glGenTextures(1, &environment_texture);
  glBindTexture(GL_TEXTURE_2D, environment_texture);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, img2.width(), img2.height(), 0, GL_BGRA,
               GL_UNSIGNED_BYTE, img2.convertToFormat(QImage::Format_ARGB32).bits());

  glGenTextures(1, &slicing_texture);
  glBindTexture(GL_TEXTURE_1D, slicing_texture);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  static const unsigned char slicing_img[] = { 0b11111111, 0b00011100 };
  glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB, 2, 0, GL_RGB, GL_UNSIGNED_BYTE_3_3_2, &slicing_img);
}

void MyViewer::draw() {
  if (model_type == ModelType::BEZIER_SURFACE && show_control_points)
    drawControlNet();

  glPolygonMode(GL_FRONT_AND_BACK, !show_solid && show_wireframe ? GL_LINE : GL_FILL);
  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1, 1);

  if (show_solid || show_wireframe) {
    if (visualization == Visualization::PLAIN)
      glColor3d(1.0, 1.0, 1.0);
    else if (visualization == Visualization::ISOPHOTES) {
      glBindTexture(GL_TEXTURE_2D, current_isophote_texture);
      glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
      glEnable(GL_TEXTURE_2D);
      glTexGeni(GL_S, GL_TEXTURE_GEN_MODE, GL_SPHERE_MAP);
      glTexGeni(GL_T, GL_TEXTURE_GEN_MODE, GL_SPHERE_MAP);
      glEnable(GL_TEXTURE_GEN_S);
      glEnable(GL_TEXTURE_GEN_T);
    } else if (visualization == Visualization::SLICING) {
      glBindTexture(GL_TEXTURE_1D, slicing_texture);
      glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
      glEnable(GL_TEXTURE_1D);
    }
    for (auto f : mesh.faces()) {
      glBegin(GL_POLYGON);
      for (auto v : mesh.fv_range(f)) {
        if (visualization == Visualization::MEAN)
          glColor3dv(meanMapColor(mesh.data(v).mean));
        else if (visualization == Visualization::SLICING)
          glTexCoord1d(mesh.point(v) | slicing_dir * slicing_scaling);
        glNormal3dv(mesh.normal(v).data());
        glVertex3dv(mesh.point(v).data());
      }
      glEnd();
    }
    if (visualization == Visualization::ISOPHOTES) {
      glDisable(GL_TEXTURE_GEN_S);
      glDisable(GL_TEXTURE_GEN_T);
      glDisable(GL_TEXTURE_2D);
      glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    } else if (visualization == Visualization::SLICING) {
      glDisable(GL_TEXTURE_1D);
    }
  }

  if (show_solid && show_wireframe) {
    glPolygonMode(GL_FRONT, GL_LINE);
    glColor3d(0.0, 0.0, 0.0);
    glDisable(GL_LIGHTING);
    for (auto f : mesh.faces()) {
      glBegin(GL_POLYGON);
      for (auto v : mesh.fv_range(f))
        glVertex3dv(mesh.point(v).data());
      glEnd();
    }
    glEnable(GL_LIGHTING);
  }

  if (axes.shown)
    drawAxes();
}

void MyViewer::drawControlNet() const {
  glDisable(GL_LIGHTING);
  glLineWidth(3.0);
  glColor3d(0.3, 0.3, 1.0);
  size_t m = degree[1] + 1;
  for (size_t k = 0; k < 2; ++k)
    for (size_t i = 0; i <= degree[k]; ++i) {
      glBegin(GL_LINE_STRIP);
      for (size_t j = 0; j <= degree[1-k]; ++j) {
        size_t const index = k ? j * m + i : i * m + j;
        const auto &p = control_points[index];
        glVertex3dv(p);
      }
      glEnd();
    }
  glLineWidth(1.0);
  glPointSize(8.0);
  glColor3d(1.0, 0.0, 1.0);
  glBegin(GL_POINTS);
  for (const auto &p : control_points)
    glVertex3dv(p);
  glEnd();
  glPointSize(1.0);
  glEnable(GL_LIGHTING);
}

void MyViewer::drawAxes() const {
  const Vec &p = axes.position;
  glColor3d(1.0, 0.0, 0.0);
  drawArrow(p, p + Vec(axes.size, 0.0, 0.0), axes.size / 50.0);
  glColor3d(0.0, 1.0, 0.0);
  drawArrow(p, p + Vec(0.0, axes.size, 0.0), axes.size / 50.0);
  glColor3d(0.0, 0.0, 1.0);
  drawArrow(p, p + Vec(0.0, 0.0, axes.size), axes.size / 50.0);
  glEnd();
}

void MyViewer::drawWithNames() {
  if (axes.shown)
    return drawAxesWithNames();

  switch (model_type) {
  case ModelType::NONE: break;
  case ModelType::MESH:
    if (!show_wireframe)
      return;
    for (auto v : mesh.vertices()) {
      glPushName(v.idx());
      glRasterPos3dv(mesh.point(v).data());
      glPopName();
    }
    break;
  case ModelType::BEZIER_SURFACE:
    if (!show_control_points)
      return;
    for (size_t i = 0, ie = control_points.size(); i < ie; ++i) {
      Vec const &p = control_points[i];
      glPushName(i);
      glRasterPos3fv(p);
      glPopName();
    }
    break;
  }
}

void MyViewer::drawAxesWithNames() const {
  const Vec &p = axes.position;
  glPushName(0);
  drawArrow(p, p + Vec(axes.size, 0.0, 0.0), axes.size / 50.0);
  glPopName();
  glPushName(1);
  drawArrow(p, p + Vec(0.0, axes.size, 0.0), axes.size / 50.0);
  glPopName();
  glPushName(2);
  drawArrow(p, p + Vec(0.0, 0.0, axes.size), axes.size / 50.0);
  glPopName();
}

void MyViewer::postSelection(const QPoint &p) {
  int sel = selectedName();
  if (sel == -1) {
    axes.shown = false;
    return;
  }

  if (axes.shown) {
    axes.selected_axis = sel;
    bool found;
    axes.grabbed_pos = camera()->pointUnderPixel(p, found);
    axes.original_pos = axes.position;
    if (!found)
      axes.shown = false;
    return;
  }

  selected_vertex = sel;
  if (model_type == ModelType::MESH)
    axes.position = Vec(mesh.point(MyMesh::VertexHandle(sel)).data());
  if (model_type == ModelType::BEZIER_SURFACE)
    axes.position = control_points[sel];
  double depth = camera()->projectedCoordinatesOf(axes.position)[2];
  Vec q1 = camera()->unprojectedCoordinatesOf(Vec(0.0, 0.0, depth));
  Vec q2 = camera()->unprojectedCoordinatesOf(Vec(width(), height(), depth));
  axes.size = (q1 - q2).norm() / 10.0;
  axes.shown = true;
  axes.selected_axis = -1;
}

void MyViewer::keyPressEvent(QKeyEvent *e) {
  if (e->modifiers() == Qt::NoModifier)
    switch (e->key()) {
    case Qt::Key_R:
      if (model_type == ModelType::MESH)
        openMesh(last_filename, false);
      else if (model_type == ModelType::BEZIER_SURFACE)
        openBezier(last_filename, false);
      update();
      break;
    case Qt::Key_O:
      if (camera()->type() == qglviewer::Camera::PERSPECTIVE)
        camera()->setType(qglviewer::Camera::ORTHOGRAPHIC);
      else
        camera()->setType(qglviewer::Camera::PERSPECTIVE);
      update();
      break;
    case Qt::Key_P:
      visualization = Visualization::PLAIN;
      update();
      break;
    case Qt::Key_M:
      visualization = Visualization::MEAN;
      update();
      break;
    case Qt::Key_L:
      visualization = Visualization::SLICING;
      update();
      break;
    case Qt::Key_I:
      visualization = Visualization::ISOPHOTES;
      current_isophote_texture = isophote_texture;
      update();
      break;
    case Qt::Key_E:
      visualization = Visualization::ISOPHOTES;
      current_isophote_texture = environment_texture;
      update();
      break;
    case Qt::Key_C:
      show_control_points = !show_control_points;
      update();
      break;
    case Qt::Key_S:
      show_solid = !show_solid;
      update();
      break;
    case Qt::Key_W:
      show_wireframe = !show_wireframe;
      update();
      break;
    case Qt::Key_F:
      fairMesh();
      update();
      break;
    case Qt::Key_U:
      elevateU();
      updateMesh();
      update();
      break;
    case Qt::Key_V:
      elevateV();
      updateMesh();
      update();
      break;
    case Qt::Key_T:
      trigonometric_basis = true;
      updateMesh();
      update();
      break;
    case Qt::Key_B:
      trigonometric_basis = false;
      updateMesh();
      update();
      break;
    default:
      QGLViewer::keyPressEvent(e);
    }
  else if (e->modifiers() == Qt::KeypadModifier)
    switch (e->key()) {
    case Qt::Key_Plus:
      slicing_scaling *= 2;
      update();
      break;
    case Qt::Key_Minus:
      slicing_scaling /= 2;
      update();
      break;
    case Qt::Key_Asterisk:
      slicing_dir = Vector(static_cast<double *>(camera()->viewDirection()));
      update();
      break;
    } else
    QGLViewer::keyPressEvent(e);
}

Vec MyViewer::intersectLines(const Vec &ap, const Vec &ad, const Vec &bp, const Vec &bd) {
  // always returns a point on the (ap, ad) line
  double a = ad * ad, b = ad * bd, c = bd * bd;
  double d = ad * (ap - bp), e = bd * (ap - bp);
  if (a * c - b * b < 1.0e-7)
    return ap;
  double s = (b * e - c * d) / (a * c - b * b);
  return ap + s * ad;
}

static void bernstein(size_t n, double u, std::vector<double> &coeff) {
  coeff.clear(); coeff.reserve(n + 1);
  coeff.push_back(1.0);
  double u1 = 1.0 - u;
  for (size_t j = 1; j <= n; ++j) {
    double saved = 0.0;
    for (size_t k = 0; k < j; ++k) {
      double tmp = coeff[k];
      coeff[k] = saved + tmp * u1;
      saved = tmp * u;
    }
    coeff.push_back(saved);
  }
}

static void bernstein(size_t n, double u, size_t derivatives,
                      std::vector<std::vector<double>> &coeffs) {
  // Assumes derivatives <= n
  std::vector<double> coeff;
  bernstein(n, u, coeff);
  coeffs.clear();
  coeffs.push_back(coeff);
  if (derivatives == 0)
    return;

  std::vector<std::vector<double>> rec;
  bernstein(n - 1, u, derivatives - 1, rec);
  for (size_t i = 1; i <= derivatives; ++i) {
    const auto &last = rec[i-1];
    std::vector<double> tmp;
    tmp.push_back(n * -last[0]);
    for (size_t j = 1; j < n; ++j)
      tmp.push_back(n * (last[j-1] - last[j]));
    tmp.push_back(n * last[n-1]);
    coeffs.push_back(tmp);
  }
}

Vec MyViewer::evaluate(double u, double v, size_t derivatives,
                       std::vector<std::vector<Vec>> &der) {
  size_t n = degree[0], m = degree[1];
  std::vector<std::vector<double>> coeff_u, coeff_v;
  if (trigonometric_basis && n % 2 == 1 && m % 2 == 1) {
    // TODO
  } else {
    bernstein(n, u, derivatives, coeff_u);
    bernstein(m, v, derivatives, coeff_v);
  }
  der.resize(derivatives + 1);
  for (size_t i = 0; i <= derivatives; ++i)
    der[i] = std::vector<Vec>(derivatives + 1, Vec(0, 0, 0));
  for (size_t i = 0; i <= derivatives; ++i)
    for (size_t j = 0; j <= derivatives; ++j)
      for (size_t k = 0, index = 0; k <= n; ++k)
        for (size_t l = 0; l <= m; ++l, ++index)
          der[i][j] += control_points[index] * coeff_u[i][k] * coeff_v[j][l];
  return der[0][0];
}

Vec MyViewer::evaluate(double u, double v) {
  static std::vector<std::vector<Vec>> dummy;
  return evaluate(u, v, 0, dummy);
}

void MyViewer::generateMesh(size_t resolution) {
  mesh.clear();
  std::vector<MyMesh::VertexHandle> handles, tri;

  for (size_t i = 0; i < resolution; ++i) {
    double u = (double)i / (double)(resolution - 1);
    for (size_t j = 0; j < resolution; ++j) {
      double v = (double)j / (double)(resolution - 1);
      auto p = evaluate(u, v);
      handles.push_back(mesh.add_vertex(Vector(static_cast<double *>(p))));
      mesh.data(handles.back()).u = u;
      mesh.data(handles.back()).v = v;
    }
  }
  for (size_t i = 0; i < resolution - 1; ++i)
    for (size_t j = 0; j < resolution - 1; ++j) {
      tri.clear();
      tri.push_back(handles[i * resolution + j]);
      tri.push_back(handles[i * resolution + j + 1]);
      tri.push_back(handles[(i + 1) * resolution + j]);
      mesh.add_face(tri);
      tri.clear();
      tri.push_back(handles[(i + 1) * resolution + j]);
      tri.push_back(handles[i * resolution + j + 1]);
      tri.push_back(handles[(i + 1) * resolution + j + 1]);
      mesh.add_face(tri);
    }
}

void MyViewer::elevateU() {
  std::vector<Vec> tmp;
  for (size_t j = 0; j <= degree[1]; ++j)
    tmp.push_back(control_points[j]);
  for (size_t i = 1, index = degree[1] + 1; i <= degree[0]; ++i) {
    double ratio = (double)i / (degree[0] + 1);
    for (size_t j = 0; j <= degree[1]; ++j, ++index)
      tmp.push_back(control_points[index-degree[1]-1] * ratio +
                    control_points[index] * (1 - ratio));
  }
  for (size_t j = 0, base = degree[0] * (degree[1] + 1); j <= degree[1]; ++j)
    tmp.push_back(control_points[base+j]);
  control_points = tmp;
  degree[0]++;
}

void MyViewer::elevateV() {
  std::vector<Vec> tmp;
  for (size_t i = 0, index = 0; i <= degree[0]; ++i) {
    tmp.push_back(control_points[index++]);
    for (size_t j = 1; j <= degree[1]; ++j, ++index) {
      double ratio = (double)j / (degree[1] + 1);
      tmp.push_back(control_points[index-1] * ratio +
                    control_points[index] * (1 - ratio));
    }
    tmp.push_back(control_points[index-1]);
  }
  control_points = tmp;
  degree[1]++;
}

void MyViewer::mouseMoveEvent(QMouseEvent *e) {
  if (!axes.shown ||
      (axes.selected_axis < 0 && !(e->modifiers() & Qt::ControlModifier)) ||
      !(e->modifiers() & (Qt::ShiftModifier | Qt::ControlModifier)) ||
      !(e->buttons() & Qt::LeftButton))
    return QGLViewer::mouseMoveEvent(e);

  if (e->modifiers() & Qt::ControlModifier) {
    // move in screen plane
    double depth = camera()->projectedCoordinatesOf(axes.position)[2];
    axes.position = camera()->unprojectedCoordinatesOf(Vec(e->pos().x(), e->pos().y(), depth));
  } else {
    Vec from, dir, axis(axes.selected_axis == 0, axes.selected_axis == 1, axes.selected_axis == 2);
    camera()->convertClickToLine(e->pos(), from, dir);
    auto p = intersectLines(axes.grabbed_pos, axis, from, dir);
    float d = (p - axes.grabbed_pos) * axis;
    axes.position[axes.selected_axis] = axes.original_pos[axes.selected_axis] + d;
  }

  if (model_type == ModelType::MESH)
    mesh.set_point(MyMesh::VertexHandle(selected_vertex),
                   Vector(static_cast<double *>(axes.position)));
  if (model_type == ModelType::BEZIER_SURFACE)
    control_points[selected_vertex] = axes.position;
  updateMesh();
  update();
}

QString MyViewer::helpString() const {
  QString text("<h2>Sample Framework</h2>"
               "<p>This is a minimal framework for 3D mesh manipulation, which can be "
               "extended and used as a base for various projects, for example "
               "prototypes for fairing algorithms, or even displaying/modifying "
               "parametric surfaces, etc.</p>"
               "<p>The following hotkeys are available:</p>"
               "<ul>"
               "<li>&nbsp;R: Reload model</li>"
               "<li>&nbsp;O: Toggle orthographic projection</li>"
               "<li>&nbsp;P: Set plain map (no coloring)</li>"
               "<li>&nbsp;M: Set mean curvature map</li>"
               "<li>&nbsp;L: Set slicing map<ul>"
               "<li>&nbsp;+: Increase slicing density</li>"
               "<li>&nbsp;-: Decrease slicing density</li>"
               "<li>&nbsp;*: Set slicing direction to view</li></ul></li>"
               "<li>&nbsp;I: Set isophote line map</li>"
               "<li>&nbsp;E: Set environment texture</li>"
               "<li>&nbsp;C: Toggle control polygon visualization</li>"
               "<li>&nbsp;S: Toggle solid (filled polygon) visualization</li>"
               "<li>&nbsp;W: Toggle wireframe visualization</li>"
               "<li>&nbsp;F: Fair mesh</li>"
               "<li>&nbsp;U: Elevate U degree (Bézier surface)</li>"
               "<li>&nbsp;V: Elevate V degree (Bézier surface)</li>"
               "<li>&nbsp;T: Change to trigonometric basis (only odd degrees)</li>"
               "<li>&nbsp;B: Change to Bernstein basis</li>"
               "</ul>"
               "<p>There is also a simple selection and movement interface, enabled "
               "only when the wireframe/controlnet is displayed: a mesh vertex can be selected "
               "by shift-clicking, and it can be moved by shift-dragging one of the "
               "displayed axes. Pressing ctrl enables movement in the screen plane.</p>"
               "<p>Note that libQGLViewer is furnished with a lot of useful features, "
               "such as storing/loading view positions, or saving screenshots. "
               "OpenMesh also has a nice collection of tools for mesh manipulation: "
               "decimation, subdivision, smoothing, etc. These can provide "
               "good comparisons to the methods you implement.</p>"
               "<p>This software can be used as a sample GUI base for handling "
               "parametric or procedural surfaces, as well. The power of "
               "Qt and libQGLViewer makes it easy to set up a prototype application. "
               "Feel free to modify and explore!</p>"
               "<p align=\"right\">Peter Salvi</p>");
  return text;
}
