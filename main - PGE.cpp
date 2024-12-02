// 3D graphics - part 4
//
// youtube video:  https://www.youtube.com/watch?v=nBzCS-Y0FcY&t=82s

// This is a continuation of 3D graphics - part 3.cbp

#include <fstream>
#include <strstream>
#include <algorithm>

//using namespace std;

#define OLC_PGE_APPLICATION
#include "olcPixelGameEngine.h"

#include       "vec3d.h"
#include      "mat4x4.h"
#include "graphics_3D.h"

struct mesh {
    std::vector<triangle> tris;

    bool LoadFromObjectFile( std::string sFilename ) {
        bool resFlag;

        std::ifstream f(sFilename);
        if (!f.is_open()) {
            resFlag = false;
        } else {
            resFlag = true;

            std::vector<vec3d> verts; // local cache of vertices

            while (!f.eof()) {
                char line[128];  // ASSUMPTION - that lines will not exceed this much characters
                f.getline( line, 128 );
                // turn the line into a stringstream for convenient processing
                std::strstream s;
                s << line;

                char junk;
                if (line[0] == 'v') {
                    vec3d v;                          // read the vertex and put it
                    s >> junk >> v.x >> v.y >> v.z;   // in the local cache of vertices
                    verts.push_back(v);
                }
                if (line[0] == 'f') {  // build triangle mesh using cache of vertices
                    int f[3];
                    s >> junk >> f[0] >> f[1] >> f[2];
                    // face indexes in object file start from 1, and not from 0. Hence - 1
                    tris.push_back({ verts[f[0] - 1], verts[f[1] - 1], verts[f[2] - 1] });
                }
                // ignore all lines not starting with a 'v' or 'f' :)
            }
        }
        return resFlag;
    }
};

class olcEngine3D : public olc::PixelGameEngine {
public:
    olcEngine3D() {
        sAppName = "3D demo - part 4";
    }

private:
    mesh meshCube;    // a std::vector of triangles, each containing three vec3d's
                      // (which in turn each contain three floats for x, y, z)
    mesh camCube;     // this cube is wrapped around camera position for testing collisions of meshes

    float fTheta = 0.0f;     // Spins world transform

    camera cCam;      // just one camera (for now)

    olc::Sprite *sprTex1;

// ==============================/   Rendering code    /==============================

// NOTE: the camera is implemented in a separate class. However everything w.r.t. rendering should be
// done in this olcConsoleGameEngine class because I want the camera to be independently implemented.

    float *pDepthBuffer = nullptr;

    void InitDepthBuffer() {
        // Depth buffer: every pixel on the screen has an associated floating point depth value
        pDepthBuffer = new float[ScreenWidth() * ScreenHeight()];
    }

    void ClearDepthBuffer() {
        for (int i = 0; i < ScreenWidth() * ScreenHeight(); i++)
            pDepthBuffer[i] = 0.0f;
    }

    // Kind of clear screen for the viewport associated with the camera
    void ClearCameraViewPort( camera &cam, bool viewPortBorder ) {
        int x1 = cam.nViewPortX1;
        int y1 = cam.nViewPortY1;
        int x2 = cam.nViewPortX2;
        int y2 = cam.nViewPortY2;

        FillRect( x1, y1, x2 - x1, y2 - y1, olc::BLACK );
        if (viewPortBorder)
            for (int i = 0; i < 3; i++) {
                DrawLine( x1 - i - 1, y1 - i - 1, x2 + i,     y1 - i - 1, olc::YELLOW );
                DrawLine( x1 - i - 1, y1 - i - 1, x1 - i - 1, y2 + i,     olc::YELLOW );
                DrawLine( x2 + i,     y2 + i,     x1 - i - 1, y2 + i,     olc::YELLOW );
                DrawLine( x2 + i,     y2 + i,     x2 + i,     y1 - i - 1, olc::YELLOW );
            }
    }

    // utility function used in DisplayCameraSettings to display radians in degrees in [0, 360)
    float radiansToDegrees( float angleInRadians ) {
        float result = angleInRadians * (180.0f / 3.14159265f);
        while (result <    0.0f) result += 360.0f;
        while (result >= 360.0f) result -= 360.0f;
        return result;
    }

    // Puts the camera setting in text on the screen.
    void DisplayCameraSettings( camera &cam, int x, int y ) {
        // test code to display camera position, orientation angles and viewport settings
        DrawString( x, y +  0, "Camera position (x):  " + std::to_string( cam.vPosition.x ),                     olc::WHITE  );
        DrawString( x, y +  1, "Camera position (y):  " + std::to_string( cam.vPosition.y ),                     olc::WHITE  );
        DrawString( x, y +  2, "Camera position (z):  " + std::to_string( cam.vPosition.z ),                     olc::WHITE  );
        DrawString( x, y +  4, "Camera angle (pitch): " + std::to_string( radiansToDegrees( cam.fCameraPitch )), olc::YELLOW );
        DrawString( x, y +  5, "Camera angle (yaw)  : " + std::to_string( radiansToDegrees( cam.fCameraYaw   )), olc::YELLOW );
        DrawString( x, y +  6, "Camera angle (roll) : " + std::to_string( radiansToDegrees( cam.fCameraRoll  )), olc::YELLOW );
        DrawString( x, y +  8, "Viewport - X1       : " + std::to_string( cam.nViewPortX1     ),                 olc::CYAN   );
        DrawString( x, y +  9, "Viewport - Y1       : " + std::to_string( cam.nViewPortY1     ),                 olc::CYAN   );
        DrawString( x, y + 10, "Viewport - X2       : " + std::to_string( cam.nViewPortX2     ),                 olc::CYAN   );
        DrawString( x, y + 11, "Viewport - Y2       : " + std::to_string( cam.nViewPortY2     ),                 olc::CYAN   );
        DrawString( x, y + 12, "Viewport - width    : " + std::to_string( cam.nViewPortWidth  ),                 olc::CYAN   );
        DrawString( x, y + 13, "Viewport - height   : " + std::to_string( cam.nViewPortHeight ),                 olc::CYAN   );
    }

    void RenderTriangles( std::vector<triangle> &trisToRender ) {

        switch (glbRenderMode) {
            case RM_TEXTURED:
                for (auto &t : trisToRender) {
                    TexturedTriangle( t.p[0].x, t.p[0].y, t.t[0].u, t.t[0].v, t.t[0].w,
                                      t.p[1].x, t.p[1].y, t.t[1].u, t.t[1].v, t.t[1].w,
                                      t.p[2].x, t.p[2].y, t.t[2].u, t.t[2].v, t.t[2].w, sprTex1 );
                }
                break;
            case RM_TEXTURED_PLUS:
                for (auto &t : trisToRender) {
                    TexturedTriangle( t.p[0].x, t.p[0].y, t.t[0].u, t.t[0].v, t.t[0].w,
                                      t.p[1].x, t.p[1].y, t.t[1].u, t.t[1].v, t.t[1].w,
                                      t.p[2].x, t.p[2].y, t.t[2].u, t.t[2].v, t.t[2].w, sprTex1 );
                    DrawTriangle(     t.p[0].x, t.p[0].y,
                                      t.p[1].x, t.p[1].y,
                                      t.p[2].x, t.p[2].y, RM_FRAMECOL_PGE );
                }
                break;
            case RM_GREYFILLED:
                for (auto &t : trisToRender ) {
                    FillTriangle( t.p[0].x, t.p[0].y, t.p[1].x, t.p[1].y, t.p[2].x, t.p[2].y, t.col );
                }
                break;
            case RM_GREYFILLED_PLUS:
                for (auto &t : trisToRender ) {
                    FillTriangle( t.p[0].x, t.p[0].y, t.p[1].x, t.p[1].y, t.p[2].x, t.p[2].y, t.col );
                    DrawTriangle( t.p[0].x, t.p[0].y, t.p[1].x, t.p[1].y, t.p[2].x, t.p[2].y, RM_FRAMECOL_PGE );
                }
                break;
            case RM_WIREFRAME:
                for (auto &t : trisToRender ) {
                    DrawTriangle( t.p[0].x, t.p[0].y, t.p[1].x, t.p[1].y, t.p[2].x, t.p[2].y, RM_FRAMECOL_PGE );
                }
                break;
        }

    }

	void TexturedTriangle( int x1, int y1, float u1, float v1, float w1,
						   int x2, int y2, float u2, float v2, float w2,
						   int x3, int y3, float u3, float v3, float w3, olc::Sprite *tex) {

        // sort input parameters in low to high y value
		if (y2 < y1) {
			std::swap(y1, y2); std::swap(x1, x2);
			std::swap(u1, u2); std::swap(v1, v2); std::swap(w1, w2);

		}
		if (y3 < y1) {
			std::swap(y1, y3); std::swap(x1, x3);
			std::swap(u1, u3); std::swap(v1, v3); std::swap(w1, w3);

		}
		if (y3 < y2) {
			std::swap(y2, y3); std::swap(x2, x3);
			std::swap(u2, u3); std::swap(v2, v3); std::swap(w2, w3);
		}
		// pixels are in integer domain, but texture coords remain floats
		// calculate convenient gradient values
		int   dy1 = y2 - y1;
		int   dx1 = x2 - x1;
		float dv1 = v2 - v1;
		float du1 = u2 - u1;
		float dw1 = w2 - w1;

		int   dy2 = y3 - y1;
		int   dx2 = x3 - x1;
		float dv2 = v3 - v1;
		float du2 = u3 - u1;
		float dw2 = w3 - w1;

		float tex_u, tex_v, tex_w;

		// NOTE: day_step and dby_step are always 1 (because pixel space
        // do the same voor u and v step values
        float dax_step = 0, dbx_step = 0,
			  du1_step = 0, dv1_step = 0,
		      du2_step = 0, dv2_step = 0,
			  dw1_step = 0, dw2_step = 0;

		if (dy1) dax_step = dx1 / (float)abs(dy1);
		if (dy2) dbx_step = dx2 / (float)abs(dy2);

		if (dy1) du1_step = du1 / (float)abs(dy1);
		if (dy1) dv1_step = dv1 / (float)abs(dy1);
		if (dy1) dw1_step = dw1 / (float)abs(dy1);

		if (dy2) du2_step = du2 / (float)abs(dy2);
		if (dy2) dv2_step = dv2 / (float)abs(dy2);
		if (dy2) dw2_step = dw2 / (float)abs(dy2);

		// start to scanline fill the triangle
		if (dy1) {
			for (int i = y1; i <= y2; i++) {
				int ax = x1 + (float)(i - y1) * dax_step;
				int bx = x1 + (float)(i - y1) * dbx_step;

				float tex_su = u1 + (float)(i - y1) * du1_step;
				float tex_sv = v1 + (float)(i - y1) * dv1_step;
				float tex_sw = w1 + (float)(i - y1) * dw1_step;

				float tex_eu = u1 + (float)(i - y1) * du2_step;
				float tex_ev = v1 + (float)(i - y1) * dv2_step;
				float tex_ew = w1 + (float)(i - y1) * dw2_step;

				if (ax > bx) {
					std::swap(ax, bx);
					std::swap(tex_su, tex_eu);
					std::swap(tex_sv, tex_ev);
					std::swap(tex_sw, tex_ew);
				}

				tex_u = tex_su;
				tex_v = tex_sv;
				tex_w = tex_sw;

				float tstep = 1.0f / ((float)(bx - ax));
				float t = 0.0f;

				for (int j = ax; j < bx; j++) {
					tex_u = (1.0f - t) * tex_su + t * tex_eu;
					tex_v = (1.0f - t) * tex_sv + t * tex_ev;
					tex_w = (1.0f - t) * tex_sw + t * tex_ew;
					if (tex_w > pDepthBuffer[i*ScreenWidth() + j]) {
                        Draw(j, i, tex->Sample(tex_u / tex_w, tex_v / tex_w));
                        pDepthBuffer[i * ScreenWidth() + j] = tex_w;
					}
					t += tstep;
				}
			}
		}

		// second half of the triangle
		dy1 = y3 - y2;
		dx1 = x3 - x2;
		dv1 = v3 - v2;
		du1 = u3 - u2;
		dw1 = w3 - w2;

		if (dy1) dax_step = dx1 / (float)abs(dy1);
		if (dy2) dbx_step = dx2 / (float)abs(dy2);

		du1_step = 0, dv1_step = 0;
		if (dy1) du1_step = du1 / (float)abs(dy1);
		if (dy1) dv1_step = dv1 / (float)abs(dy1);
		if (dy1) dw1_step = dw1 / (float)abs(dy1);

		if (dy1) {
			for (int i = y2; i <= y3; i++) {
				int ax = x2 + (float)(i - y2) * dax_step;
				int bx = x1 + (float)(i - y1) * dbx_step;

				float tex_su = u2 + (float)(i - y2) * du1_step;
				float tex_sv = v2 + (float)(i - y2) * dv1_step;
				float tex_sw = w2 + (float)(i - y2) * dw1_step;

				float tex_eu = u1 + (float)(i - y1) * du2_step;
				float tex_ev = v1 + (float)(i - y1) * dv2_step;
				float tex_ew = w1 + (float)(i - y1) * dw2_step;

				if (ax > bx) {
					std::swap(ax, bx);
					std::swap(tex_su, tex_eu);
					std::swap(tex_sv, tex_ev);
					std::swap(tex_sw, tex_ew);
				}

				tex_u = tex_su;
				tex_v = tex_sv;
				tex_w = tex_sw;

				float tstep = 1.0f / ((float)(bx - ax));
				float t = 0.0f;

				for (int j = ax; j < bx; j++) {
					tex_u = (1.0f - t) * tex_su + t * tex_eu;
					tex_v = (1.0f - t) * tex_sv + t * tex_ev;
					tex_w = (1.0f - t) * tex_sw + t * tex_ew;
					if (tex_w > pDepthBuffer[i*ScreenWidth() + j]) {
                        Draw(j, i, tex->Sample(tex_u / tex_w, tex_v / tex_w));
                        pDepthBuffer[i * ScreenWidth() + j] = tex_w;
					}
					t += tstep;
				}
			}
		}
	}

// ==============================/   End of rendering code    /==============================

public:
    bool OnUserCreate() override {
        // create the depth buffer
        InitDepthBuffer();

        // populate the meshCube using an object file.
//        meshCube.LoadFromObjectFile("mountains.obj");

        // Initialize the unit cube, including texturing coordinates
		meshCube.tris = {
            { 0.0f, 0.0f, 0.0f, 1.0f,   0.0f, 1.0f, 0.0f, 1.0f,   1.0f, 1.0f, 0.0f, 1.0f,   0.0f, 1.0f, 1.0f,   0.0f, 0.0f, 1.0f,   1.0f, 0.0f, 1.0f, },    // SOUTH
            { 0.0f, 0.0f, 0.0f, 1.0f,   1.0f, 1.0f, 0.0f, 1.0f,   1.0f, 0.0f, 0.0f, 1.0f,   0.0f, 1.0f, 1.0f,   1.0f, 0.0f, 1.0f,   1.0f, 1.0f, 1.0f, },
            { 1.0f, 0.0f, 0.0f, 1.0f,   1.0f, 1.0f, 0.0f, 1.0f,   1.0f, 1.0f, 1.0f, 1.0f,   0.0f, 1.0f, 1.0f,   0.0f, 0.0f, 1.0f,   1.0f, 0.0f, 1.0f, },    // EAST
            { 1.0f, 0.0f, 0.0f, 1.0f,   1.0f, 1.0f, 1.0f, 1.0f,   1.0f, 0.0f, 1.0f, 1.0f,   0.0f, 1.0f, 1.0f,   1.0f, 0.0f, 1.0f,   1.0f, 1.0f, 1.0f, },
            { 1.0f, 0.0f, 1.0f, 1.0f,   1.0f, 1.0f, 1.0f, 1.0f,   0.0f, 1.0f, 1.0f, 1.0f,   0.0f, 1.0f, 1.0f,   0.0f, 0.0f, 1.0f,   1.0f, 0.0f, 1.0f, },    // NORTH
            { 1.0f, 0.0f, 1.0f, 1.0f,   0.0f, 1.0f, 1.0f, 1.0f,   0.0f, 0.0f, 1.0f, 1.0f,   0.0f, 1.0f, 1.0f,   1.0f, 0.0f, 1.0f,   1.0f, 1.0f, 1.0f, },
            { 0.0f, 0.0f, 1.0f, 1.0f,   0.0f, 1.0f, 1.0f, 1.0f,   0.0f, 1.0f, 0.0f, 1.0f,   0.0f, 1.0f, 1.0f,   0.0f, 0.0f, 1.0f,   1.0f, 0.0f, 1.0f, },    // WEST
            { 0.0f, 0.0f, 1.0f, 1.0f,   0.0f, 1.0f, 0.0f, 1.0f,   0.0f, 0.0f, 0.0f, 1.0f,   0.0f, 1.0f, 1.0f,   1.0f, 0.0f, 1.0f,   1.0f, 1.0f, 1.0f, },
            { 0.0f, 1.0f, 0.0f, 1.0f,   0.0f, 1.0f, 1.0f, 1.0f,   1.0f, 1.0f, 1.0f, 1.0f,   0.0f, 1.0f, 1.0f,   0.0f, 0.0f, 1.0f,   1.0f, 0.0f, 1.0f, },    // TOP
            { 0.0f, 1.0f, 0.0f, 1.0f,   1.0f, 1.0f, 1.0f, 1.0f,   1.0f, 1.0f, 0.0f, 1.0f,   0.0f, 1.0f, 1.0f,   1.0f, 0.0f, 1.0f,   1.0f, 1.0f, 1.0f, },
            { 1.0f, 0.0f, 1.0f, 1.0f,   0.0f, 0.0f, 1.0f, 1.0f,   0.0f, 0.0f, 0.0f, 1.0f,   0.0f, 1.0f, 1.0f,   0.0f, 0.0f, 1.0f,   1.0f, 0.0f, 1.0f, },    // BOTTOM
            { 1.0f, 0.0f, 1.0f, 1.0f,   0.0f, 0.0f, 0.0f, 1.0f,   1.0f, 0.0f, 0.0f, 1.0f,   0.0f, 1.0f, 1.0f,   1.0f, 0.0f, 1.0f,   1.0f, 1.0f, 1.0f, },
		};

		// make a wrapping cube for the camera position to test collision detection methods
		camCube = meshCube;

//		sprTex1 = new olcSprite( "ceiling_tiled.org.spr" );
		sprTex1 = new olc::Sprite( "Jario.png" );

        // Create cameras defining the associated viewports
        float sw = (float)ScreenWidth();
        float sh = (float)ScreenHeight();
        cCam.InitCamera( this, "the camera", (int)(sw * 0.05f), (int)(sh * 0.05f), (int)(sw * 0.95f), (int)(sh * 0.95f), 60.0f, 0.1f, 1000.0f);

        FillRect( 0, 0, ScreenWidth(), ScreenHeight(), olc::DARK_GREEN );

        fTheta = 0.0f;

        return true;
    }

// ======================/ COLLISION DETECTION MESH vs MESH  /======================

// see: https://docs.godotengine.org/en/stable/tutorials/math/vectors_advanced.html#collision-detection-in-3d

    // this function analyses collision between two meshes (sets of triangles)
    // NOTE: the sets of triangles must be transformed into world space to do
    // meaningful collision detection
    // NOTE: the function appears to be working correctly only for convex objects
    bool AreObjectsColliding( std::vector<triangle> &vObjA, std::vector<triangle> &vObjB ) {

        bool bOverlapping = true;

        // In this first loop all vertices of mesh B are checked against the triangles of mesh A
        for (auto triA : vObjA) {   // iterate over all planes (triangles) of vObjA
            bool bAllOut = true;
            // Calculate the normal vector from the current triangle and normalize it
            vec3d normal = Vector_GetNormal( triA.p[0], triA.p[1], triA.p[2] );
            // iterate over all vertices of vObjB
            for (auto triB: vObjB) {
                for (int i = 0; bAllOut && i < 3; i++) {
                    // calculate distance of point of triB to plane of triA
                    if (Vector_Distance( normal, triA.p[0], triB.p[i] ) < 0.0f) bAllOut = false;
                }
                if (!bAllOut) break;
            }
            if (bAllOut) {
                bOverlapping = false;
                break;
            }
        }

        // Now the other way around: all vertices of mesh A are checked against the triangles of mesh B
        if (bOverlapping) {
            for (auto triB : vObjB) {   // iterate over all planes (triangles) of vObjB
                bool bAllOut = true;
                // Calculate the normal vector from the current triangle and normalize it
                vec3d normal = Vector_GetNormal( triB.p[0], triB.p[1], triB.p[2] );
                // iterate over all vertices of vObjA
                for (auto triA: vObjA) {
                    for (int i = 0; bAllOut && i < 3; i++) {
                        // calculate distance of point of triB to plane of triA
                        if (Vector_Distance( normal, triB.p[0], triA.p[i] ) < 0.0f) bAllOut = false;
                    }
                    if (!bAllOut) break;
                }
                if (bAllOut) {
                    bOverlapping = false;
                    break;
                }
            }
        }

        if (bOverlapping) {
            for (auto triA : vObjA) {
                for (int i = 0; i < 3; i++) {   // iterate all edges of vObjA in ea
                    vec3d ea = Vector_Sub( triA.p[i], triA.p[(i + 1) % 3] );

                    for (auto triB : vObjB) {
                        for (int j = 0; j < 3; j++) {   // iterate all edges of vObjB in eb
                            vec3d eb = Vector_Sub( triB.p[j], triB.p[(j + 1) % 3] );

                            vec3d n = Vector_CrossProduct( ea, eb );
                            if (Vector_Length( n ) != 0.0f) {

                                float max_A = -INFINITY;
                                float min_A =  INFINITY;
                                for (auto tA : vObjA) {
                                    for (int a = 0; a < 3; a++) {
                                        float d = Vector_DotProduct( n, tA.p[a] );
                                        max_A = std::max( max_A, d );
                                        min_A = std::min( min_A, d );
                                    }
                                }

                                float max_B = -INFINITY;
                                float min_B =  INFINITY;
                                for (auto tB : vObjB) {
                                    for (int b = 0; b < 3; b++) {
                                        float d = Vector_DotProduct( n, tB.p[b] );
                                        max_B = std::max( max_B, d );
                                        min_B = std::min( min_B, d );
                                    }
                                }

                                if (min_A > max_B || min_B > max_A) {
                                    bOverlapping = false;
                                }
                            }
                        }
                        if (!bOverlapping)
                            break;
                    }
                }
                if (!bOverlapping)
                    break;
            }
        }

        return bOverlapping;
    }

    bool OnUserUpdate( float fElapsedTime ) {

// ====================/ User input / ====================

        if (GetKey( olc::Key::K1 ).bHeld) glbRenderMode = RM_TEXTURED;
        if (GetKey( olc::Key::K2 ).bHeld) glbRenderMode = RM_TEXTURED_PLUS;
        if (GetKey( olc::Key::K3 ).bHeld) glbRenderMode = RM_GREYFILLED;
        if (GetKey( olc::Key::K4 ).bHeld) glbRenderMode = RM_GREYFILLED_PLUS;
        if (GetKey( olc::Key::K5 ).bHeld) glbRenderMode = RM_WIREFRAME;

// First get the angles updated because some of the other movements depend on it.

        // adapt yaw angle (around x axis)
        if (GetKey( olc::Key::Z ).bHeld) cCam.fCameraPitch -= (GetKey( olc::Key::SHIFT ).bHeld ? 0.8f : 0.2f) * fElapsedTime;
        if (GetKey( olc::Key::Q ).bHeld) cCam.fCameraPitch += (GetKey( olc::Key::SHIFT ).bHeld ? 0.8f : 0.2f) * fElapsedTime;// adapt yaw angle (around y axis)
        if (GetKey( olc::Key::A ).bHeld) cCam.fCameraYaw   -= (GetKey( olc::Key::SHIFT ).bHeld ? 0.8f : 0.2f) * fElapsedTime;
        if (GetKey( olc::Key::D ).bHeld) cCam.fCameraYaw   += (GetKey( olc::Key::SHIFT ).bHeld ? 0.8f : 0.2f) * fElapsedTime;// adapt yaw angle (around z axis)
        if (GetKey( olc::Key::E ).bHeld) cCam.fCameraRoll  -= (GetKey( olc::Key::SHIFT ).bHeld ? 0.8f : 0.2f) * fElapsedTime;
        if (GetKey( olc::Key::R ).bHeld) cCam.fCameraRoll  += (GetKey( olc::Key::SHIFT ).bHeld ? 0.8f : 0.2f) * fElapsedTime;

// ====================/ Camera updates / ====================

        cCam.RecalculateCamera();

// ====================/ More user input / ====================

        // strafe forward / backwards in the direction of the vLookDir vector of the camera.
        vec3d vLookTmp = Vector_Mul( cCam.vLookDir, (GetKey( olc::Key::SHIFT ).bHeld ? 4.0f : 1.0f) * fElapsedTime );
        if (GetKey( olc::Key::W ).bHeld) cCam.vPosition = Vector_Add( cCam.vPosition, vLookTmp );
        if (GetKey( olc::Key::S ).bHeld) cCam.vPosition = Vector_Sub( cCam.vPosition, vLookTmp );
        // strafe up / down in the direction of the vUp vector of the camera.
        vec3d vUpTmp = Vector_Mul( cCam.vUp, (GetKey( olc::Key::SHIFT ).bHeld ? 4.0f : 1.0f) * fElapsedTime );
        if (GetKey( olc::Key::UP   ).bHeld) cCam.vPosition = Vector_Add( cCam.vPosition, vUpTmp );
        if (GetKey( olc::Key::DOWN ).bHeld) cCam.vPosition = Vector_Sub( cCam.vPosition, vUpTmp );
        // strafe left / right perpendicular to the direction the camera is looking, using the vRight vector.
        vec3d vRightTmp = Vector_Mul( cCam.vRight, (GetKey( olc::Key::SHIFT ).bHeld ? 4.0f : 1.0f) * fElapsedTime );
        if (GetKey( olc::Key::LEFT  ).bHeld) cCam.vPosition = Vector_Sub( cCam.vPosition, vRightTmp );
        if (GetKey( olc::Key::RIGHT ).bHeld) cCam.vPosition = Vector_Add( cCam.vPosition, vRightTmp );

// ====================/ Start 3D graphics calculations / ====================

// shut the rotating off for testing
//        fTheta += 1.0f * fElapsedTime;

		// Store triangles for rastering later
		std::vector<triangle> vecWorldCubeTris,
                         vecWorldCamTris,
                         vecTrianglesToRaster,
                         vecTrianglesToRender;

        // Create the world matrix for the transform to world coordinates, using scaling parameters (first three)
        // the rotation angles around x, y and z axis and translation distances in x, y and z direction (last three parameters).
        mat4x4 matWorldCube = Matrix_MakeTransformComplete( 1.0f, 1.0f, 1.0f, fTheta, 0.0f, fTheta * 0.5f, 0.0f, 0.0f, 6.0f );

        // Prepare triangles for rendering
        for (auto triOriginal : meshCube.tris) {
            triangle triTransformed;

            // For each point in the triangle, transform it using the world matrix.
            // The effect is the rotation around z- and x-axes, and the translation
			// (= offset into the screen - to make sure the object is actually in front of viewer)

			// Note that this function is part of the camera object, but needs only
            // be called once - EVEN IN THE CASE OF MULTIPLE CAMERAS!
            cCam.Tri_WorldTransform( triOriginal, matWorldCube, triTransformed );
            vecWorldCubeTris.push_back( triTransformed );
        }

        mat4x4 matWorldCam = Matrix_MakeTransformComplete( 1.0f, 1.0f, 1.0f,   // no up/down scaling
                                                           0.0f, 0.0f, 0.0f,   // no rotation angles necessary
                                                           cCam.vPosition.x - 0.5f,   // translate so that origin of camera
                                                           cCam.vPosition.y - 0.5f,   // coincides with center of camCube
                                                           cCam.vPosition.z - 0.5f );
        for (auto triOriginal : camCube.tris) {
            triangle triTransformed;

            cCam.Tri_WorldTransform( triOriginal, matWorldCam, triTransformed );
            vecWorldCamTris.push_back( triTransformed );
        }
        bool bCollDetected = AreObjectsColliding( vecWorldCubeTris, vecWorldCamTris );

        for (auto triTransformed : vecWorldCubeTris) {

            // Do the culling, the view and project transform per camera. The output is added
            // to the vector that is passed as parameter.
            // NOTE: clipping against near plane is done in this function.
//            vec3d vLight = { cCam.vLookDir.x, cCam.vLookDir.y, - cCam.vLookDir.z };
            cCam.CullViewAndProjectTriangle( triTransformed, vecTrianglesToRaster );
        }
        // do the clipping against the borders of the viewport and produce a list to render
        cCam.RasterizeTriangles( vecTrianglesToRaster, vecTrianglesToRender );

        if (bCollDetected) {
            FillRect( 0, 0, ScreenWidth(), ScreenHeight(), olc::DARK_MAGENTA );
        } else {
            FillRect( 0, 0, ScreenWidth(), ScreenHeight(), olc::DARK_GREEN );
        }

        // Clear depth buffer
        ClearDepthBuffer();
        // Clear viewports
        ClearCameraViewPort( cCam, true );
        // finally render the results
        RenderTriangles( vecTrianglesToRender );

        return true;
    }
};

// Vary the pixel size between 2 and 8. The High def constants are the max screen width / height for
// pixel size == 1.
// NOTE: a pixel size < 4 is deteriorated performance, and a pixel size > 6 is degraded imagery
#define PIXEL_SIZE      1

#define HIGH_DEF_X   1024
#define HIGH_DEF_Y    600

int main()
{
    olcEngine3D demo;
    if (demo.Construct( HIGH_DEF_X / PIXEL_SIZE, HIGH_DEF_Y / PIXEL_SIZE, PIXEL_SIZE, PIXEL_SIZE ))
        demo.Start();
    else
        std::cout << "ERROR: couldn't start, because demo.ConstructConsole() call returned false" << std::endl;

    return 0;
}
