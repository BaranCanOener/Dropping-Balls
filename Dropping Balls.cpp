#define OLC_PGE_APPLICATION
#define _USE_MATH_DEFINES
#include "olcPixelGameEngine.h"
#include <string>
#include <time.h>
#include <thread>
#include <chrono>
#include <math.h>
#include <complex>

enum class ObstacleType { OBS_BALL, OBS_LINE , OBS_UNDEFINED};

enum class SettingParameter { SET_RESTITUTION, SET_GRAVITY, SET_DRAG, SET_SECONDS_TO_PREDICT, SET_SIM_PRECISION };

enum class AppState { STATE_MENU, STATE_PAUSE, STATE_SANDBOX, STATE_GAME};

enum class GameState {GAMESTATE_AIM, GAMESTATE_SHOOTING, GAMESTATE_NEWOBSTACLES, GAMESTATE_END};

const float DEFAULT_RADIUS = 0.1f;

const int DEFAULT_OBS_HEALTH = 10;
const int DEFAULT_BALL_HEALTH = 7;

//Coordinate system - screen edges
const float LEFT_EDGE = 0.0f;
const float RIGHT_EDGE = 10.0f;
const float TOP_EDGE = 0.0f;
const float BOTTOM_EDGE = -12.5f;

const float BALL_DENSITY = 1.0f;

struct Ball
{
	int health = DEFAULT_BALL_HEALTH;
	int ID = 0;
	bool invincible = true;
	float px = 0, py = 0;
	float vx = 0, vy = 0;
	float ax = 0, ay = 0;
	float radius = DEFAULT_BALL_HEALTH;
	float mass = 0;
	olc::Pixel colour = olc::WHITE;
};

class Obstacle
{
public:
	
	int health = DEFAULT_OBS_HEALTH;
	int ID = 0;
	bool invincible = true;
	olc::Pixel colour = olc::RED;
	ObstacleType type = ObstacleType::OBS_UNDEFINED;
};

class BallObstacle : public Obstacle
{
public:
	float px, py;
	float radius;

	BallObstacle(float px, float py, float radius, int health, bool invincible)
	{
		this->radius = radius;
		this->px = px;
		this->py = py;
		this->health = health;
		this->invincible = invincible;
		type = ObstacleType::OBS_BALL;
	}
};

class LineObstacle : public Obstacle
{
public:
	float sx, sy;
	float ex, ey;
	float radius;
	bool isFrame = false;

	LineObstacle(float sx, float sy, float ex, float ey, float radius, int health, bool invincible)
	{
		this->sx = sx;
		this->sy = sy;
		this->ex = ex;
		this->ey = ey;
		this->radius = radius;
		this->health = health;
		this->invincible = invincible;
		type = ObstacleType::OBS_LINE;
	}

};

class PhysicsSystem
{
public:
	//Physics parameters
	float gravity = 0.0f; //meters per second
	float restitutionFactor = 1.0f; //controls elasticity of collisions
	float drag = 0.0f;
	bool timeStop = false;
	bool invincibleObstacles = false;
	bool ballVballCollisions = true;

	int shapesSpawned = 0;

	//Coordinates of spawning balls
	float spawn_x = 5.0f;
	float spawn_y = 0.0f;

	std::vector<Ball> balls;
	std::vector<Obstacle*> obstacles;

	//Keep track of the score (i.e. the number of hits of balls against obstacles)
	int score = 0;

	//Simulation Parameters - how far into the future a balls flight path is simulated
	float simulation_seconds_to_predict = 3.0f;
	float simulation_timedelta = 0.005f;
	int simulation_timesteps = (int)(simulation_seconds_to_predict / simulation_timedelta);

	//Deep copies the entire data from the argument into the class; i.e. the systems parameters and objects are fully replicated
	void copyFrom(PhysicsSystem system)
	{
		gravity = system.gravity;
		restitutionFactor = system.restitutionFactor;
		drag = system.drag;
		timeStop = system.timeStop;
		invincibleObstacles = system.invincibleObstacles;
		ballVballCollisions = system.ballVballCollisions;
		shapesSpawned = system.shapesSpawned;
		spawn_x = system.spawn_x;
		spawn_y = system.spawn_y;
		score = system.score;
		simulation_seconds_to_predict = system.simulation_seconds_to_predict;
		simulation_timedelta = system.simulation_timedelta;
		simulation_timesteps = system.simulation_timesteps;
		balls = system.balls;

		while (!obstacles.empty())
		{
			delete obstacles.back();
			obstacles.pop_back();
		}


		for (int i = 0; i < system.obstacles.size(); i++)
		{
			if (system.obstacles.at(i)->type == ObstacleType::OBS_LINE)
			{
				LineObstacle* orig_obstacle = static_cast<LineObstacle*>(system.obstacles.at(i));
				LineObstacle* cpy_obstacle = new LineObstacle(orig_obstacle->sx, orig_obstacle->sy, orig_obstacle->ex, orig_obstacle->ey, orig_obstacle->radius, orig_obstacle->health, orig_obstacle->invincible);
				cpy_obstacle->ID = orig_obstacle->ID;
				obstacles.push_back(cpy_obstacle);
			}
			if (system.obstacles.at(i)->type == ObstacleType::OBS_BALL)
			{
				BallObstacle* orig_obstacle = static_cast<BallObstacle*>(system.obstacles.at(i));
				BallObstacle* cpy_obstacle = new BallObstacle(orig_obstacle->px, orig_obstacle->py, orig_obstacle->radius, orig_obstacle->health, orig_obstacle->invincible);
				cpy_obstacle->ID = orig_obstacle->ID;
				obstacles.push_back(cpy_obstacle);
			}
		}
	}

	void cleanup()
	{
		//Cleanup routines
		for (unsigned int i = 0; i < balls.size(); i++)
		{
			//Remove balls that are destroyed, below the bottom, left, or right boundary, or far enough above the top boundary
			if ((balls.at(i).health <= 0) || (balls.at(i).py < BOTTOM_EDGE) || (balls.at(i).px > RIGHT_EDGE) || (balls.at(i).px < LEFT_EDGE) || (balls.at(i).py > TOP_EDGE + (TOP_EDGE - BOTTOM_EDGE)))
				balls.erase(balls.begin() + i);
		}
		for (unsigned int i = 0; i < obstacles.size(); i++)
		{
			if (obstacles.at(i)->health <= 0)
			{
				delete obstacles.at(i);
				obstacles.erase(obstacles.begin() + i);
			}

		}
	}

	void SpawnBall(float px, float py, float vx, float vy, float ax, float ay, float radius, int health, bool invincible, olc::Pixel colour)
	{
		shapesSpawned++;
		Ball ball;
		ball.ax = ax;
		ball.ay = ay;
		ball.px = px;
		ball.py = py;
		ball.vx = vx;
		ball.vy = vy;
		ball.radius = radius;
		ball.health = health;
		ball.invincible = invincible;
		ball.mass = M_PI * radius * radius * BALL_DENSITY;
		ball.ID = shapesSpawned;
		ball.colour = colour;
		balls.push_back(ball);
	}

	BallObstacle* SpawnObstacle_Ball(float px, float py, float radius, int health, bool invincible)
	{
		BallObstacle* obstacle = new BallObstacle(px, py, radius, health, invincible);
		shapesSpawned++;
		obstacle->ID = shapesSpawned;
		obstacles.push_back(obstacle);
		return obstacle;
	}

	LineObstacle* SpawnObstacle_Line(float sx, float sy, float ex, float ey, float radius, int health, bool invincible)
	{
		LineObstacle* obstacle = new LineObstacle(sx, sy, ex, ey, radius, health, invincible);
		shapesSpawned++;
		obstacle->ID = shapesSpawned;
		obstacles.push_back(obstacle);
		return obstacle;
	}

	bool ballsOverlap(float px1, float py1, float r1, float px2, float py2, float r2)
	{
		return (px1 - px2) * (px1 - px2) + (py1 - py2) * (py1 - py2) < (r1 + r2) * (r1 + r2);
	}

	bool staticDisplaceBallBall(Ball* ball1, Ball* ball2)
	{
		float distance = sqrtf((ball1->px - ball2->px) * (ball1->px - ball2->px) + (ball1->py - ball2->py) * (ball1->py - ball2->py));
		float overlap = 0.5f * (distance - ball1->radius - ball2->radius);

		if (overlap <= 0) {
			ball1->px -= overlap * (ball1->px - ball2->px) / distance;
			ball1->py -= overlap * (ball1->py - ball2->py) / distance;

			ball2->px += overlap * (ball1->px - ball2->px) / distance;
			ball2->py += overlap * (ball1->py - ball2->py) / distance;
			return true;
		}
		else
			return false;
	}

	//Handles the collision between a nonstationary ball and a stationary obstacle.
	//Returns true if ball collides with obstacle, statically resolves the collision,
	//and calls dynamicResolveBallBall with a dummy ball to control the velocity change of ball 
	//prevLocation: Will be used for more precise collision detection
	bool dynamicResolveBallObstacle(Ball* ball, Obstacle* obstacle, Ball* prevLocation)
	{
		bool collision = false;
		Ball* collisionBall = new Ball();
		collisionBall->vx = -ball->vx;
		collisionBall->vy = -ball->vy;
		collisionBall->mass = ball->mass;
		if (obstacle->type == ObstacleType::OBS_BALL)
		{
			float distance = sqrtf((ball->px - static_cast<BallObstacle*>(obstacle)->px) * (ball->px - static_cast<BallObstacle*>(obstacle)->px) + (ball->py - static_cast<BallObstacle*>(obstacle)->py) * (ball->py - static_cast<BallObstacle*>(obstacle)->py));
			float overlap = (distance - ball->radius - static_cast<BallObstacle*>(obstacle)->radius);
			if (distance <= (ball->radius + static_cast<BallObstacle*>(obstacle)->radius))
			{
				collision = true;
				ball->px -= overlap * (ball->px - static_cast<BallObstacle*>(obstacle)->px) / distance;
				ball->py -= overlap * (ball->py - static_cast<BallObstacle*>(obstacle)->py) / distance;
				collisionBall->px = static_cast<BallObstacle*>(obstacle)->px;
				collisionBall->py = static_cast<BallObstacle*>(obstacle)->py;
				dynamicResolveBallBall(ball, collisionBall);
			}

		}
		else if (obstacle->type == ObstacleType::OBS_LINE)
		{
			LineObstacle* lObs = static_cast<LineObstacle*>(obstacle);
			float lineX1 = lObs->ex - lObs->sx;
			float lineY1 = lObs->ey - lObs->sy;

			float lineX2 = ball->px - lObs->sx;
			float lineY2 = ball->py - lObs->sy;

			float lineLength = lineX1 * lineX1 + lineY1 * lineY1;

			float s = max(0, min(lineLength, (lineX1 * lineX2 + lineY1 * lineY2))) / lineLength;

			float closestPointX = lObs->sx + s * lineX1;
			float closestPointY = lObs->sy + s * lineY1;

			float distance = sqrt((ball->px - closestPointX) * (ball->px - closestPointX) + (ball->py - closestPointY) * (ball->py - closestPointY));

			if (distance <= (ball->radius + lObs->radius))
			{	
				collision = true;
				collisionBall->vx = -ball->vx;
				collisionBall->vy = -ball->vy;
				collisionBall->px = closestPointX;
				collisionBall->py = closestPointY;
				collisionBall->radius = lObs->radius;
				float overlap = (distance - ball->radius - collisionBall->radius);
				ball->px -= overlap * (ball->px - collisionBall->px) / distance;
				ball->py -= overlap * (ball->py - collisionBall->py) / distance;
				dynamicResolveBallBall(ball, collisionBall);
			}
		}
		delete collisionBall;
		return collision;
	}

	void dynamicResolveBallBall(Ball* ball1, Ball* ball2)
	{
		float distance = sqrtf((ball1->px - ball2->px) * (ball1->px - ball2->px) + (ball1->py - ball2->py) * (ball1->py - ball2->py));

		//normal vector
		float nx = (ball1->px - ball2->px) / distance;
		float ny = (ball1->py - ball2->py) / distance;

		//tangent vector
		float tx = -ny;
		float ty = nx;

		//For the tangential response: dot product between ball's velocity and tangent vector
		float dpTan = ball1->vx * tx + ball1->vy * ty;
		float dpTan2 = ball2->vx * tx + ball2->vy * ty;

		//For the normal response: dot product between the ball's velocity and normal vector
		float dpNorm = ball1->vx * nx + ball1->vy * ny;
		float dpNorm2 = ball2->vx * nx + ball2->vy * ny;

		//Conservation of momentum in 1D
		float m1 = (dpNorm * (ball1->mass - ball2->mass) + 2.0f * ball2->mass * dpNorm2) / (ball1->mass + ball2->mass);
		float m2 = (dpNorm2 * (ball2->mass - ball1->mass) + 2.0f * ball1->mass * dpNorm) / (ball1->mass + ball2->mass);

		ball1->vx = (tx * dpTan + nx * m1) * restitutionFactor;
		ball1->vy = (ty * dpTan + ny * m1) * restitutionFactor;
		ball2->vx = (tx * dpTan2 + nx * m2) * restitutionFactor;
		ball2->vy = (ty * dpTan2 + ny * m2) * restitutionFactor;
	}


	//fixedResolution = false: invokes a single simulation step, i.e. the resolution depends on the frames per second
	//fixedResolution = true: invokes simulation steps as dictated by simulation_timedelta until fElapsedTime is reached.
	//						  In this case, the simulation will overshoot by (ceil(fElapsedTime/simulation_timedelta) * simulation_timedelta - fElapsedTime,
	//						  i.e. in case of high fps counts, the movements per frame will be too fast.
	//						  One fix was to make the thread sleep in order to cap the frame counter, but a more precise solution should be implemented.
	float physicsHandling(float fElapsedTime, bool fixedResolution)
	{
		float simTimestep;

		float simTotalElapsedTime = 0.0f;
		if (!fixedResolution)
			simTimestep = fElapsedTime;
		else
			simTimestep = simulation_timedelta;
		//Collision detection between balls v balls; and balls v obstacles
		std::vector<std::tuple<Ball*, Ball*>> collidingPairs;
		while (simTotalElapsedTime < fElapsedTime)
		{
			simTotalElapsedTime += simTimestep;

			for (auto& ball : balls)
			{
				Ball prevLocation = ball;

				//Controls movement of balls
				ball.ay = -gravity - ball.vy * drag;
				ball.ax = -ball.vx * drag;
				ball.vx += ball.ax * simTimestep;
				ball.vy += ball.ay * simTimestep;
				ball.px += ball.vx * simTimestep;
				ball.py += ball.vy * simTimestep;
				if (ballVballCollisions)
				{
					for (auto& ball2 : balls) {
						if (ball.ID != ball2.ID)
						{
							//Statically resolve collision - displace balls so they don't overlap
							if (staticDisplaceBallBall(&ball, &ball2))
							{
								//Dynamically resolve collision, i.e. alter velocity vectors
								dynamicResolveBallBall(&ball, &ball2); 
								if (!ball.invincible)
									ball.health--;
								if (!ball2.invincible)
									ball2.health--;
							}
						}
					}
				}

				//handle collision between balls and obstacles
				for (auto& obstacle : obstacles)
				{
					//Resolve collision, i.e. (a) displace ball [obstacle is static] and (b) alter the balls velocity vector
					if (dynamicResolveBallObstacle(&ball, obstacle, &prevLocation))
					{
						if ((!obstacle->invincible) && (!invincibleObstacles))
						{
							obstacle->health -= 1;
							score++;
						}
						if (!ball.invincible)
							ball.health--;
					}
				}

			}
		}
		return simTotalElapsedTime;
	}
};

class AIPlayer
{
public:
	PhysicsSystem physics;

	float secondsPerSimulatedShot = 8.0f;

	//Returns the score for a simulated shot at coordinate (x,y) in the physics-System
	int simulateShot(float x, float y, int ballsToSpawn, float timeBetweenBalls)
	{
		float seconds = 0.0f;
		float timeSinceSpawn = 0.0f;
		physics.score = 0;
		while (seconds < secondsPerSimulatedShot)
		{
			physics.physicsHandling(physics.simulation_timedelta, false);

			//Cleanup routines
			physics.cleanup();

			timeSinceSpawn += physics.simulation_timedelta;
			if ((ballsToSpawn > 0) && (timeSinceSpawn > timeBetweenBalls))
			{
				timeSinceSpawn = 0;
				ballsToSpawn--;
				physics.SpawnBall(physics.spawn_x, physics.spawn_y, (x - physics.spawn_x) * 2, (y - physics.spawn_y) * 2, 0, 0, DEFAULT_RADIUS, DEFAULT_BALL_HEALTH, false, olc::Pixel(rand() % 256, rand() % 256, rand() % 256));
			}

			seconds += physics.simulation_timedelta;
		}
		return physics.score;
	}
};

class BallsDrop : public olc::PixelGameEngine
{
private:

	AppState appState = AppState::STATE_MENU;
	GameState gameState = GameState::GAMESTATE_NEWOBSTACLES;
	PhysicsSystem physics;
	AIPlayer AI;

	//For conversion between world- to screen coordinates and vice versa
	float pixelsPerUnit = 0.0f;

	//Runtime
	float simulationTime = 0.0f;

	//Round-Parameters
	int round = 0;
	const int spawnAfterRounds = 1; //Controls after how many rounds a new set of obstacles spawn
	const float minMove = 1.0f; //Controls the minimum distance of an up-move when new obstacles spawn
	const float maxMove = 2.0f; //Controls the maximum distance of an up-move when new obstacles spawn
	const float moveUpInSeconds = 1.0f; //Controls the time it takes for an up-move when new obstacles spawn

	//Variables to control the shooting of new balls
	int ballsPerShot = 10; //Total amount of balls a player shoots
	int ballsToSpawn = 0;
	float timeBetweenBalls = 0.01f;
	float timeSinceSpawn = 0.0f;
	std::tuple<float, float> targetCoord;
	
	//Controls the current up-move
	float Y_move = 0.0f; 
	float seconds_moved = 0.0f;

	//Controls the amount of obstacles that spawn each round
	const int minObstacleSpawn = 3;
	const int maxObstacleSpawn = 5;

	//Drawing of debug menu
	bool debugMenu = true;
	
	Obstacle* selectedObstacle;
	ObstacleType spawnMode = ObstacleType::OBS_BALL;
	SettingParameter setting = SettingParameter::SET_GRAVITY;

	void spawnInputHandling()
	{
		std::tuple<float, float> wCoord;
		std::tuple<int, int> sCoord;

		//Spawn Balls
		if (GetMouse(0).bReleased)
		{
			std::get<0>(sCoord) = GetMouseX();
			std::get<1>(sCoord) = GetMouseY();
			wCoord = S_to_W(sCoord);
			physics.SpawnBall(physics.spawn_x, physics.spawn_y, (std::get<0>(wCoord) - physics.spawn_x) * 2, (std::get<1>(wCoord) - physics.spawn_y) * 2, 0, 0, DEFAULT_RADIUS, DEFAULT_BALL_HEALTH, true, olc::Pixel(rand() % 256, rand() % 256, rand() % 256));
		}

		//Spawn Obstacles
		if (spawnMode == ObstacleType::OBS_BALL)
		{
			if (GetMouse(1).bReleased)
			{
				std::get<0>(sCoord) = GetMouseX();
				std::get<1>(sCoord) = GetMouseY();
				wCoord = S_to_W(sCoord);
				physics.SpawnObstacle_Ball(std::get<0>(wCoord), std::get<1>(wCoord), DEFAULT_RADIUS * 2, DEFAULT_OBS_HEALTH, false);
			}
		}
		else if (spawnMode == ObstacleType::OBS_LINE)
		{
			if (GetMouse(1).bPressed)
			{
				std::get<0>(sCoord) = GetMouseX();
				std::get<1>(sCoord) = GetMouseY();
				wCoord = S_to_W(sCoord);
				selectedObstacle = physics.SpawnObstacle_Line(std::get<0>(wCoord), std::get<1>(wCoord), std::get<0>(wCoord), std::get<1>(wCoord), 0.001f, DEFAULT_OBS_HEALTH, false);
			}
			if ((GetMouse(1).bHeld) || (GetMouse(1).bReleased))
			{
				std::get<0>(sCoord) = GetMouseX();
				std::get<1>(sCoord) = GetMouseY();
				wCoord = S_to_W(sCoord);
				static_cast<LineObstacle*>(selectedObstacle)->ex = std::get<0>(wCoord);
				static_cast<LineObstacle*>(selectedObstacle)->ey = std::get<1>(wCoord);
			}
		}
	}

	void spawnFrame()
	{
		physics.SpawnObstacle_Line(LEFT_EDGE, TOP_EDGE, LEFT_EDGE, BOTTOM_EDGE, 0.01f, DEFAULT_OBS_HEALTH, true);
		static_cast<LineObstacle*>(physics.obstacles.back())->isFrame = true;
		physics.SpawnObstacle_Line(RIGHT_EDGE, TOP_EDGE, RIGHT_EDGE, BOTTOM_EDGE, 0.01f, DEFAULT_OBS_HEALTH, true);
		static_cast<LineObstacle*>(physics.obstacles.back())->isFrame = true;
		physics.SpawnObstacle_Line(LEFT_EDGE, BOTTOM_EDGE, RIGHT_EDGE, BOTTOM_EDGE, 0.01f, DEFAULT_OBS_HEALTH, true);
		static_cast<LineObstacle*>(physics.obstacles.back())->isFrame = true;
	}

	float random_float(float min, float max)
	{
		return min + (static_cast <float> (rand()) / static_cast <float> (RAND_MAX))*(max - min);
	}

	void spawnMainMenuObjects()
	{
		//Head
		physics.SpawnObstacle_Line((float)(LEFT_EDGE + (RIGHT_EDGE - LEFT_EDGE) / 5), (float)(BOTTOM_EDGE + (TOP_EDGE - BOTTOM_EDGE) / 2), (float)(LEFT_EDGE + (RIGHT_EDGE - LEFT_EDGE) * 4 / 5), (float)(BOTTOM_EDGE + (TOP_EDGE - BOTTOM_EDGE) / 2), 0.001f, DEFAULT_OBS_HEALTH, true);
		physics.SpawnObstacle_Line((float)(LEFT_EDGE + (RIGHT_EDGE - LEFT_EDGE) / 5), (float)(BOTTOM_EDGE + (TOP_EDGE - BOTTOM_EDGE) / 8), (float)(LEFT_EDGE + (RIGHT_EDGE - LEFT_EDGE) * 4 / 5), (float)(BOTTOM_EDGE + (TOP_EDGE - BOTTOM_EDGE) / 8), 0.001f, DEFAULT_OBS_HEALTH, true);
		physics.SpawnObstacle_Line((float)(LEFT_EDGE + (RIGHT_EDGE - LEFT_EDGE) / 5), (float)(BOTTOM_EDGE + (TOP_EDGE - BOTTOM_EDGE) / 2), (float)(LEFT_EDGE + (RIGHT_EDGE - LEFT_EDGE) / 5), (float)(BOTTOM_EDGE + (TOP_EDGE - BOTTOM_EDGE) / 8), 0.001f, DEFAULT_OBS_HEALTH, true);
		physics.SpawnObstacle_Line((float)(LEFT_EDGE + (RIGHT_EDGE - LEFT_EDGE) * 4 / 5), (float)(BOTTOM_EDGE + (TOP_EDGE - BOTTOM_EDGE) / 2), LEFT_EDGE + (RIGHT_EDGE - LEFT_EDGE) * 4 / 5, (float)(BOTTOM_EDGE + (TOP_EDGE - BOTTOM_EDGE) / 8), 0.001f, DEFAULT_OBS_HEALTH, true);
		
		//Left Eye
		physics.SpawnObstacle_Line((float)(LEFT_EDGE + (RIGHT_EDGE - LEFT_EDGE) * 1.5 / 5), (float)(BOTTOM_EDGE + (TOP_EDGE - BOTTOM_EDGE) / 2.25), (float)(LEFT_EDGE + (RIGHT_EDGE - LEFT_EDGE) * 2.0 / 5), (float)(BOTTOM_EDGE + (TOP_EDGE - BOTTOM_EDGE) / 2.25), 0.001f, DEFAULT_OBS_HEALTH, true);
		physics.SpawnObstacle_Line((float)(LEFT_EDGE + (RIGHT_EDGE - LEFT_EDGE) * 1.5 / 5), (float)(BOTTOM_EDGE + (TOP_EDGE - BOTTOM_EDGE) / 2.7), (float)(LEFT_EDGE + (RIGHT_EDGE - LEFT_EDGE) * 2.0 / 5), (float)(BOTTOM_EDGE + (TOP_EDGE - BOTTOM_EDGE) / 2.7), 0.001f, DEFAULT_OBS_HEALTH, true);
		physics.SpawnObstacle_Line((float)(LEFT_EDGE + (RIGHT_EDGE - LEFT_EDGE) * 1.5 / 5), (float)(BOTTOM_EDGE + (TOP_EDGE - BOTTOM_EDGE) / 2.7), (float)(LEFT_EDGE + (RIGHT_EDGE - LEFT_EDGE) * 1.5 / 5), (float)(BOTTOM_EDGE + (TOP_EDGE - BOTTOM_EDGE) / 2.25), 0.001f, DEFAULT_OBS_HEALTH, true);
		physics.SpawnObstacle_Line((float)(LEFT_EDGE + (RIGHT_EDGE - LEFT_EDGE) * 2.0 / 5), (float)(BOTTOM_EDGE + (TOP_EDGE - BOTTOM_EDGE) / 2.7), (float)(LEFT_EDGE + (RIGHT_EDGE - LEFT_EDGE) * 2.0 / 5), (float)(BOTTOM_EDGE + (TOP_EDGE - BOTTOM_EDGE) / 2.25), 0.001f, DEFAULT_OBS_HEALTH, true);
		
		//Right Eye
		physics.SpawnObstacle_Line((float)(LEFT_EDGE + (RIGHT_EDGE - LEFT_EDGE) * 3 / 5), (float)(BOTTOM_EDGE + (TOP_EDGE - BOTTOM_EDGE) / 2.25), (float)(LEFT_EDGE + (RIGHT_EDGE - LEFT_EDGE) * 3.5 / 5), (float)(BOTTOM_EDGE + (TOP_EDGE - BOTTOM_EDGE) / 2.25), 0.001f, DEFAULT_OBS_HEALTH, true);
		physics.SpawnObstacle_Line((float)(LEFT_EDGE + (RIGHT_EDGE - LEFT_EDGE) * 3 / 5), (float)(BOTTOM_EDGE + (TOP_EDGE - BOTTOM_EDGE) / 2.7), (float)(LEFT_EDGE + (RIGHT_EDGE - LEFT_EDGE) * 3.5 / 5), (float)(BOTTOM_EDGE + (TOP_EDGE - BOTTOM_EDGE) / 2.7), 0.001f, DEFAULT_OBS_HEALTH, true);
		physics.SpawnObstacle_Line((float)(LEFT_EDGE + (RIGHT_EDGE - LEFT_EDGE) * 3 / 5), (float)(BOTTOM_EDGE + (TOP_EDGE - BOTTOM_EDGE) / 2.7), (float)(LEFT_EDGE + (RIGHT_EDGE - LEFT_EDGE) * 3.0 / 5), (float)(BOTTOM_EDGE + (TOP_EDGE - BOTTOM_EDGE) / 2.25), 0.001f, DEFAULT_OBS_HEALTH, true);
		physics.SpawnObstacle_Line((float)(LEFT_EDGE + (RIGHT_EDGE - LEFT_EDGE) * 3.5 / 5), (float)(BOTTOM_EDGE + (TOP_EDGE - BOTTOM_EDGE) / 2.7), (float)(LEFT_EDGE + (RIGHT_EDGE - LEFT_EDGE) * 3.5 / 5), (float)(BOTTOM_EDGE + (TOP_EDGE - BOTTOM_EDGE) / 2.25), 0.001f, DEFAULT_OBS_HEALTH, true);
		
		//Mouth
		physics.SpawnObstacle_Line((float)(LEFT_EDGE + (RIGHT_EDGE - LEFT_EDGE) * 1.5 / 5), (float)(BOTTOM_EDGE + (TOP_EDGE - BOTTOM_EDGE) * 2 / 8), (float)(LEFT_EDGE + (RIGHT_EDGE - LEFT_EDGE) * 3.5 / 5), (float)(BOTTOM_EDGE + (TOP_EDGE - BOTTOM_EDGE) * 2 / 8), 0.001f, DEFAULT_OBS_HEALTH, true);
		physics.SpawnObstacle_Line((float)(LEFT_EDGE + (RIGHT_EDGE - LEFT_EDGE) * 1.5 / 5), (float)(BOTTOM_EDGE + (TOP_EDGE - BOTTOM_EDGE) * 1.5 / 8), (float)(LEFT_EDGE + (RIGHT_EDGE - LEFT_EDGE) * 3.5 / 5), (float)(BOTTOM_EDGE + (TOP_EDGE - BOTTOM_EDGE) * 1.5 / 8), 0.001f, DEFAULT_OBS_HEALTH, true);
		physics.SpawnObstacle_Line((float)(LEFT_EDGE + (RIGHT_EDGE - LEFT_EDGE) * 1.5 / 5), (float)(BOTTOM_EDGE + (TOP_EDGE - BOTTOM_EDGE) * 1.5 / 8), (float)(LEFT_EDGE + (RIGHT_EDGE - LEFT_EDGE) * 1.5 / 5), (float)(BOTTOM_EDGE + (TOP_EDGE - BOTTOM_EDGE) * 2 / 8), 0.001f, DEFAULT_OBS_HEALTH, true);
		physics.SpawnObstacle_Line((float)(LEFT_EDGE + (RIGHT_EDGE - LEFT_EDGE) * 2 / 5), (float)(BOTTOM_EDGE + (TOP_EDGE - BOTTOM_EDGE) * 1.5 / 8), (float)(LEFT_EDGE + (RIGHT_EDGE - LEFT_EDGE) * 2 / 5), (float)(BOTTOM_EDGE + (TOP_EDGE - BOTTOM_EDGE) * 2 / 8), 0.001f, DEFAULT_OBS_HEALTH, true);
		physics.SpawnObstacle_Line((float)(LEFT_EDGE + (RIGHT_EDGE - LEFT_EDGE) * 2.5 / 5), (float)(BOTTOM_EDGE + (TOP_EDGE - BOTTOM_EDGE) * 1.5 / 8), (float)(LEFT_EDGE + (RIGHT_EDGE - LEFT_EDGE) * 2.5 / 5), (float)(BOTTOM_EDGE + (TOP_EDGE - BOTTOM_EDGE) * 2 / 8), 0.001f, DEFAULT_OBS_HEALTH, true);
		physics.SpawnObstacle_Line((float)(LEFT_EDGE + (RIGHT_EDGE - LEFT_EDGE) * 3 / 5), (float)(BOTTOM_EDGE + (TOP_EDGE - BOTTOM_EDGE) * 1.5 / 8), (float)(LEFT_EDGE + (RIGHT_EDGE - LEFT_EDGE) * 3 / 5), (float)(BOTTOM_EDGE + (TOP_EDGE - BOTTOM_EDGE) * 2 / 8), 0.001f, DEFAULT_OBS_HEALTH, true);
		physics.SpawnObstacle_Line((float)(LEFT_EDGE + (RIGHT_EDGE - LEFT_EDGE) * 3.5 / 5), (float)(BOTTOM_EDGE + (TOP_EDGE - BOTTOM_EDGE) * 1.5 / 8), (float)(LEFT_EDGE + (RIGHT_EDGE - LEFT_EDGE) * 3.5 / 5), (float)(BOTTOM_EDGE + (TOP_EDGE - BOTTOM_EDGE) * 2 / 8), 0.001f, DEFAULT_OBS_HEALTH, true);

		std::tuple<float, float> sCoord, wCoord;

		//Spawn 25 balls inside the robot's head
		for (int i = 0; i < 25; i++)
		{
			std::get<0>(sCoord) = (float)(rand() % ScreenWidth());
			std::get<1>(sCoord) = (float)(rand() % ScreenHeight());
			wCoord = S_to_W(sCoord);
			//Linear transformation from whole world space to compressed frame defined by the above lines
			std::get<0>(wCoord) = (float)(LEFT_EDGE + (RIGHT_EDGE - LEFT_EDGE) / 5 + std::get<0>(wCoord) * (((RIGHT_EDGE - LEFT_EDGE) * 3 / 5) / (RIGHT_EDGE - LEFT_EDGE)));
			std::get<1>(wCoord) = (float)(BOTTOM_EDGE + (TOP_EDGE - BOTTOM_EDGE) / 2 + std::get<1>(wCoord) * (((TOP_EDGE - BOTTOM_EDGE) * 3 / 8) / (TOP_EDGE - BOTTOM_EDGE)));
			physics.SpawnBall(std::get<0>(wCoord), std::get<1>(wCoord), random_float(-4.0, 4.0), random_float(-4.0, 4.0), 0, 0, DEFAULT_RADIUS, DEFAULT_BALL_HEALTH, true, olc::Pixel(rand() % 256, rand() % 256, rand() % 256));
		}
	}

	//Handle the spawning of several balls with a prescribed delay timeBetweenBalls
	void massBallSpawn(float fElapsedTime)
	{
		if (!physics.timeStop)
			timeSinceSpawn += fElapsedTime;
		if ((ballsToSpawn > 0) && (timeSinceSpawn > timeBetweenBalls))
		{
			std::tuple<float, float> wCoord;
			std::tuple<int, int> sCoord;
			timeSinceSpawn = 0;
			ballsToSpawn--;
			if (appState == AppState::STATE_SANDBOX)
				physics.SpawnBall(physics.spawn_x, physics.spawn_y, (std::get<0>(targetCoord) - physics.spawn_x) * 2, (std::get<1>(targetCoord) - physics.spawn_y) * 2, 0, 0, DEFAULT_RADIUS, DEFAULT_BALL_HEALTH, true, olc::Pixel(rand() % 256, rand() % 256, rand() % 256));
			else
				physics.SpawnBall(physics.spawn_x, physics.spawn_y, (std::get<0>(targetCoord) - physics.spawn_x) * 2, (std::get<1>(targetCoord) - physics.spawn_y) * 2, 0, 0, DEFAULT_RADIUS, DEFAULT_BALL_HEALTH, false, olc::WHITE);
		}
	}

	void game_moveObstacles(float x, float y)
	{
		for (auto& obstacle : physics.obstacles)
		{
			if (obstacle->type == ObstacleType::OBS_BALL)
			{
				static_cast<BallObstacle*>(obstacle)->px += x;
				static_cast<BallObstacle*>(obstacle)->py += y;
			}
			else if (obstacle->type == ObstacleType::OBS_LINE)
			{
				if (!(static_cast<LineObstacle*>(obstacle)->isFrame))
				{
					static_cast<LineObstacle*>(obstacle)->sx += x;
					static_cast<LineObstacle*>(obstacle)->sy += y;
					static_cast<LineObstacle*>(obstacle)->ex += x;
					static_cast<LineObstacle*>(obstacle)->ey += y;
				}
			}
		}
	}

	int game_countObstacles()
	{
		int count = 0;
		for (auto& obstacle : physics.obstacles)
		{
			if (obstacle->type == ObstacleType::OBS_LINE)
				if (!(static_cast<LineObstacle*>(obstacle)->isFrame))
					count++;
		}
		return count;
	}

	void game_spawnNewObstacles()
	{
		Y_move = random_float(minMove, maxMove);
		int obstaclesToSpawn = rand() % maxObstacleSpawn + minObstacleSpawn;
		float x_step = (RIGHT_EDGE - LEFT_EDGE) / (obstaclesToSpawn);
		float sx = LEFT_EDGE;
		float sy, ex, ey;
		for (int i = 0; i < obstaclesToSpawn; i++)
		{
			int rnd = rand() % 2;
			if (rnd == 0)
			{
				sx = x_step+i + random_float(0.0f, x_step*0.9f);
				ex = sx + x_step - random_float(0.0f, x_step*0.9f);
				sy = random_float(BOTTOM_EDGE - Y_move, BOTTOM_EDGE);
				ey = random_float(BOTTOM_EDGE - Y_move, BOTTOM_EDGE);
				physics.SpawnObstacle_Line(sx, sy, ex, ey, 0.001f, 10, false);
			}
			else
			{
				sx = x_step*i + random_float(0.0f, x_step*0.9f);
				sy = random_float(BOTTOM_EDGE - Y_move, BOTTOM_EDGE);
				physics.SpawnObstacle_Ball(sx, sy, random_float(DEFAULT_RADIUS * 4, DEFAULT_RADIUS * 8), DEFAULT_OBS_HEALTH,  false);
			}
		}
	}

	void despawnBalls()
	{
		for (int i = (int)(physics.balls.size() - 1); i >= 0; i--)
		{
			physics.balls.erase(physics.balls.begin() + i);
		}
	}

	void despawnObstacles()
	{
		for (int i = (int)(physics.obstacles.size() - 1); i >= 0; i--)
		{
			delete physics.obstacles.at(i);
			physics.obstacles.erase(physics.obstacles.begin() + i);
		}
	}

	//Handles the user input to configure the simulation
	void configurationInputHandling()
	{
		if ((GetKey(olc::Key::UP).bPressed))
		{
			if (setting == SettingParameter::SET_SECONDS_TO_PREDICT)
				setting = SettingParameter::SET_DRAG;
			else if (setting == SettingParameter::SET_DRAG)
				setting = SettingParameter::SET_RESTITUTION;
			else if (setting == SettingParameter::SET_RESTITUTION)
				setting = SettingParameter::SET_GRAVITY;
			else if (setting == SettingParameter::SET_GRAVITY)
				setting = SettingParameter::SET_SIM_PRECISION;
			else
				setting = SettingParameter::SET_SECONDS_TO_PREDICT;
		}
		else if ((GetKey(olc::Key::DOWN).bPressed))
		{
			if (setting == SettingParameter::SET_GRAVITY)
				setting = SettingParameter::SET_RESTITUTION;
			else if (setting == SettingParameter::SET_RESTITUTION)
				setting = SettingParameter::SET_DRAG;
			else if (setting == SettingParameter::SET_DRAG)
				setting = SettingParameter::SET_SECONDS_TO_PREDICT;
			else if (setting == SettingParameter::SET_SECONDS_TO_PREDICT)
				setting = SettingParameter::SET_SIM_PRECISION;
			else
				setting = SettingParameter::SET_GRAVITY;
		}
		else if (GetKey(olc::Key::LEFT).bPressed)
		{
			if (setting == SettingParameter::SET_GRAVITY)
				physics.gravity -= 1;
			else if (setting == SettingParameter::SET_RESTITUTION)
				physics.restitutionFactor -= 0.01f;
			else if (setting == SettingParameter::SET_DRAG)
				physics.drag -= 0.1f;
			else if (setting == SettingParameter::SET_SECONDS_TO_PREDICT)
				physics.simulation_seconds_to_predict -= 1;
			else if (setting == SettingParameter::SET_SIM_PRECISION)
			{
				physics.simulation_timedelta -= 0.001f;
				if (physics.simulation_timedelta <= 0)
					physics.simulation_timedelta = 0.001f;
			}
		}
		else if (GetKey(olc::Key::RIGHT).bPressed)
		{
			if (setting == SettingParameter::SET_GRAVITY)
				physics.gravity += 1;
			else if (setting == SettingParameter::SET_RESTITUTION)
				physics.restitutionFactor += 0.01f;
			else if (setting == SettingParameter::SET_DRAG)
				physics.drag += 0.1f;
			else if (setting == SettingParameter::SET_SECONDS_TO_PREDICT)
				physics.simulation_seconds_to_predict += 1;
			else if (setting == SettingParameter::SET_SIM_PRECISION)
				physics.simulation_timedelta += 0.001f;
		}
		else if (GetKey(olc::Key::NP0).bPressed)
		{
			if (setting == SettingParameter::SET_GRAVITY)
				physics.gravity = 0.0;
			else if (setting == SettingParameter::SET_RESTITUTION)
				physics.restitutionFactor = 0.0;
			else if (setting == SettingParameter::SET_DRAG)
				physics.drag = 0.0;
		}
		else if (GetKey(olc::Key::I).bPressed)
		{
			physics.invincibleObstacles = !physics.invincibleObstacles;
		}
		else if (GetKey(olc::Key::M).bPressed)
		{
			despawnObstacles();
		}
		else if (GetKey(olc::Key::N).bPressed)
		{
			despawnBalls();
		}
		else if (GetKey(olc::Key::B).bPressed)
		{
			physics.ballVballCollisions = !physics.ballVballCollisions;
		}
		else if (GetKey(olc::Key::O).bPressed)
		{
			if (spawnMode == ObstacleType::OBS_BALL)
				spawnMode = ObstacleType::OBS_LINE;
			else
				spawnMode = ObstacleType::OBS_BALL;
		}
		else if (GetKey(olc::Key::X).bPressed)
		{
			physics.simulation_timedelta -= 0.00001f;
		}
		else if (GetKey(olc::Key::T).bPressed)
		{
			physics.timeStop = !physics.timeStop;
		}
		else if (GetKey(olc::Key::S).bPressed) //Spawn Mass number of obstacles
		{
			std::tuple<float, float> wCoord;
			std::tuple<int, int> sCoord;
			for (int i = 0; i < 10; i++)
			{
				std::get<0>(sCoord) = rand() % ScreenWidth();
				std::get<1>(sCoord) = rand() % ScreenHeight();
				wCoord = S_to_W(sCoord);
				physics.SpawnObstacle_Ball(std::get<0>(wCoord), std::get<1>(wCoord), (rand() / (RAND_MAX / DEFAULT_RADIUS)) + DEFAULT_RADIUS, DEFAULT_OBS_HEALTH, false);
			}
		}
		else if (GetKey(olc::Key::F).bPressed)
		{
			ballsToSpawn += 20;
			std::tuple<int, int> sCoord;
			std::get<0>(sCoord) = GetMouseX();
			std::get<1>(sCoord) = GetMouseY();
			targetCoord = S_to_W(sCoord);
		}
		else if (GetKey(olc::Key::ESCAPE).bPressed)
		{
			returnToMenu();
		}
	}

	void returnToMenu()
	{
		despawnBalls();
		despawnObstacles();
		spawnMainMenuObjects();
		physics.gravity = 0.0f;
		physics.restitutionFactor = 1.0f;
		physics.drag = 0.0f;
		physics.ballVballCollisions = true;
		appState = AppState::STATE_MENU;
	}

	void drawConfigurationMenu()
	{
		DrawString(ScreenWidth() - 140, 5, "Balls: " + std::to_string(physics.balls.size()), olc::WHITE);
		DrawString(ScreenWidth() - 140, 15, "Obstacles: " + std::to_string(physics.obstacles.size()), olc::WHITE);
		if (setting == SettingParameter::SET_GRAVITY)
			DrawString(ScreenWidth() - 140, 25, "Gravity: " + std::to_string(physics.gravity), olc::RED);

		else
			DrawString(ScreenWidth() - 140, 25, "Gravity: " + std::to_string(physics.gravity), olc::WHITE);
		if (setting == SettingParameter::SET_RESTITUTION)
			DrawString(ScreenWidth() - 140, 35, "Restitution: " + std::to_string(physics.restitutionFactor), olc::RED);
		else
			DrawString(ScreenWidth() - 140, 35, "Restitution: " + std::to_string(physics.restitutionFactor), olc::WHITE);
		if (setting == SettingParameter::SET_DRAG)
			DrawString(ScreenWidth() - 140, 45, "Drag: " + std::to_string(physics.drag), olc::RED);
		else
			DrawString(ScreenWidth() - 140, 45, "Drag: " + std::to_string(physics.drag), olc::WHITE);
		if (setting == SettingParameter::SET_SECONDS_TO_PREDICT)
			DrawString(ScreenWidth() - 140, 55, "Forecast: " + std::to_string(physics.simulation_seconds_to_predict) + "s", olc::RED);
		else
			DrawString(ScreenWidth() - 140, 55, "Forecast: " + std::to_string(physics.simulation_seconds_to_predict) + "s", olc::WHITE);
		if (setting == SettingParameter::SET_SIM_PRECISION)
			DrawString(ScreenWidth() - 140, 65, "Precision: " + std::to_string(physics.simulation_timedelta) + "s", olc::RED);
		else
			DrawString(ScreenWidth() - 140, 65, "Precision: " + std::to_string(physics.simulation_timedelta) + "s", olc::WHITE);
		std::tuple<float, float, float> momentumEnergy = systemMomentum_and_energy();
		DrawString(ScreenWidth() - 140, 75, "Kinetic Energy: ", olc::WHITE);
		DrawString(ScreenWidth() - 120, 85, " " + std::to_string(std::get<2>(momentumEnergy)) + " J", olc::WHITE);
		DrawString(ScreenWidth() - 140, 95, "Momentum: ", olc::WHITE);
		DrawString(ScreenWidth() - 120, 105, "(" + std::to_string(std::get<0>(momentumEnergy)), olc::WHITE);
		DrawString(ScreenWidth() - 120, 115, " " + std::to_string(std::get<1>(momentumEnergy)) + ")", olc::WHITE);
		DrawString(ScreenWidth() - 140, 125, "Sim Time: " + std::to_string(simulationTime), olc::WHITE);
		DrawString(ScreenWidth() - 160, ScreenHeight() - 150, "KEYS:", olc::WHITE);
		DrawString(ScreenWidth() - 160, ScreenHeight() - 140, "Esc: Main menu,", olc::WHITE);
		DrawString(ScreenWidth() - 160, ScreenHeight() - 130, "d: Close menu,", olc::WHITE);
		DrawString(ScreenWidth() - 160, ScreenHeight() - 120, "t: Timestop,", olc::WHITE);
		DrawString(ScreenWidth() - 160, ScreenHeight() - 110, "LMouse: Spawn Ball,", olc::WHITE);
		DrawString(ScreenWidth() - 160, ScreenHeight() - 100, "f: Spawn 20 Balls,", olc::WHITE);
		DrawString(ScreenWidth() - 160, ScreenHeight() - 90, "RMouse: Spawn Obstacle,", olc::WHITE);
		DrawString(ScreenWidth() - 160, ScreenHeight() - 80, "Up/down: select,", olc::WHITE);
		DrawString(ScreenWidth() - 160, ScreenHeight() - 70, "Left/right: set;", olc::WHITE);
		DrawString(ScreenWidth() - 160, ScreenHeight() - 60, "n: del balls,", olc::WHITE);
		DrawString(ScreenWidth() - 160, ScreenHeight() - 50, "m: del obstcl;", olc::WHITE);
		if (physics.invincibleObstacles)
			DrawString(ScreenWidth() - 160, ScreenHeight() - 40, "i: invinc. obstc", olc::GREEN);
		else
			DrawString(ScreenWidth() - 160, ScreenHeight() - 40, "i: invinc. obstc", olc::WHITE);
		if (physics.ballVballCollisions)
			DrawString(ScreenWidth() - 160, ScreenHeight() - 30, "b: ball/ball coll.", olc::GREEN);
		else
			DrawString(ScreenWidth() - 160, ScreenHeight() - 30, "b: ball/ball coll.", olc::WHITE);
		DrawString(ScreenWidth() - 160, ScreenHeight() - 20, "o: spawn mode;", olc::WHITE);
		if (spawnMode == ObstacleType::OBS_BALL)
			DrawString(ScreenWidth() - 160, ScreenHeight() - 10, " => Ball (click)", olc::RED);
		else
			DrawString(ScreenWidth() - 160, ScreenHeight() - 10, " => Line (hold&drag)", olc::RED);
	}

	void drawSimulatedBallPath(float sx, float sy, float vx, float vy, bool invincibleBall, olc::Pixel colour)
	{
		std::tuple<float, float> wCoord, wCoord2;
		std::tuple<int, int> sCoord, sCoord2;
		physics.SpawnBall(sx, sy, vx, vy, 0, 0, DEFAULT_RADIUS, DEFAULT_BALL_HEALTH, invincibleBall, olc::WHITE);
		physics.shapesSpawned--;
		std::get<0>(wCoord) = physics.spawn_x;
		std::get<1>(wCoord) = physics.spawn_y;
		sCoord = W_to_S(wCoord);
		Ball* ball = &physics.balls.back();
		physics.simulation_timesteps = (int)(physics.simulation_seconds_to_predict / physics.simulation_timedelta);
		for (int i = 0; i < physics.simulation_timesteps; i++)
		{
			ball->ay = -physics.gravity - ball->vy * physics.drag;
			ball->ax = -ball->vx * physics.drag;
			ball->vx += ball->ax * physics.simulation_timedelta;
			ball->vy += ball->ay * physics.simulation_timedelta;
			ball->px += ball->vx * physics.simulation_timedelta;
			ball->py += ball->vy * physics.simulation_timedelta;
			for (auto& obstacle : physics.obstacles)
			{
				if (physics.dynamicResolveBallObstacle(ball, obstacle, ball))
				{
					if (!ball->invincible)
						ball->health--;
				}
			}
			std::get<0>(wCoord2) = ball->px;
			std::get<1>(wCoord2) = ball->py;
			sCoord2 = W_to_S(wCoord2);
			if (ball->health > 0)
				DrawLine(std::get<0>(sCoord), std::get<1>(sCoord), std::get<0>(sCoord2), std::get<1>(sCoord2), colour);
			sCoord = sCoord2;
		}
		physics.balls.pop_back();
	}

	olc::Pixel invertColour(olc::Pixel colour)
	{
		return olc::Pixel(255 - colour.r, 255 - colour.g, 255 - colour.b);
	}

	//Fades between green for health = DEFAULT_OBS_HEALTH and red for health = 1
	olc::Pixel healthToColour(int health)
	{
		return olc::Pixel((uint8_t)(255 * (float)(DEFAULT_OBS_HEALTH - health + 1) / (float)DEFAULT_OBS_HEALTH), (uint8_t)(255 * ((float)health / (float)DEFAULT_OBS_HEALTH)), 0);
	}

	void drawObjects()
	{
		std::tuple<float, float> wCoord, wCoord2;
		std::tuple<int, int> sCoord, sCoord2;
		//Draw Balls
		for (auto& ball : physics.balls)
		{
			std::get<0>(wCoord) = ball.px;
			std::get<1>(wCoord) = ball.py;
			sCoord = W_to_S(wCoord);
			FillCircle(std::get<0>(sCoord), std::get<1>(sCoord), int(ball.radius * pixelsPerUnit), ball.colour);
		}

		//Draw Obstacles
		for (auto& obstacle : physics.obstacles)
		{
			if (obstacle->health > 0)
			{
				olc::Pixel colour = healthToColour(obstacle->health);
				if (obstacle->type == ObstacleType::OBS_BALL)
				{
					std::get<0>(wCoord) = static_cast<BallObstacle*>(obstacle)->px;
					std::get<1>(wCoord) = static_cast<BallObstacle*>(obstacle)->py;
					sCoord = W_to_S(wCoord);
					if (appState == AppState::STATE_SANDBOX)
						FillCircle(std::get<0>(sCoord), std::get<1>(sCoord), int(static_cast<BallObstacle*>(obstacle)->radius * pixelsPerUnit), colour);
					else
						DrawCircle(std::get<0>(sCoord), std::get<1>(sCoord), int(static_cast<BallObstacle*>(obstacle)->radius * pixelsPerUnit), colour);
					if ((!obstacle->invincible) && (!physics.invincibleObstacles))
						DrawString(std::get<0>(sCoord) - 7, std::get<1>(sCoord) - 3, std::to_string(obstacle->health), olc::WHITE);
				}
				if (obstacle->type == ObstacleType::OBS_LINE)
				{
					LineObstacle* line = static_cast<LineObstacle*>(obstacle);
					std::get<0>(wCoord) = line->sx;
					std::get<1>(wCoord) = line->sy;
					std::get<0>(wCoord2) = line->ex;
					std::get<1>(wCoord2) = line->ey;
					sCoord = W_to_S(wCoord);
					sCoord2 = W_to_S(wCoord2);
					FillCircle(std::get<0>(sCoord), std::get<1>(sCoord), int(line->radius * pixelsPerUnit), colour);
					FillCircle(std::get<0>(sCoord2), std::get<1>(sCoord2), int(line->radius * pixelsPerUnit), colour);
					float nx = -(line->ey - line->sy);
					float ny = (line->ex - line->sx);
					float d = sqrt(nx * nx + ny * ny);
					nx /= d;
					ny /= d;
					std::get<0>(wCoord) = line->sx + nx * line->radius;
					std::get<1>(wCoord) = line->sy + ny * line->radius;
					std::get<0>(wCoord2) = line->ex + nx * line->radius;
					std::get<1>(wCoord2) = line->ey + ny * line->radius;
					sCoord = W_to_S(wCoord);
					sCoord2 = W_to_S(wCoord2);
					DrawLine(std::get<0>(sCoord), std::get<1>(sCoord), std::get<0>(sCoord2), std::get<1>(sCoord2), colour);
					std::get<0>(wCoord) = line->sx - nx * line->radius;
					std::get<1>(wCoord) = line->sy - ny * line->radius;
					std::get<0>(wCoord2) = line->ex - nx * line->radius;
					std::get<1>(wCoord2) = line->ey - ny * line->radius;
					sCoord = W_to_S(wCoord);
					sCoord2 = W_to_S(wCoord2);
					DrawLine(std::get<0>(sCoord), std::get<1>(sCoord), std::get<0>(sCoord2), std::get<1>(sCoord2), colour);
					if ((!obstacle->invincible) && (!physics.invincibleObstacles))
						DrawString(std::get<0>(sCoord) + (int)(0.5f * (std::get<0>(sCoord2) - std::get<0>(sCoord))), std::get<1>(sCoord) + (int)(0.5f * (std::get<1>(sCoord2) - std::get<1>(sCoord))), std::to_string(obstacle->health), olc::WHITE);
				}
			}
		}
	}

	std::tuple<float, float, float> systemMomentum_and_energy()
	{
		std::tuple<float, float, float> momentum_and_energy;
		std::get<0>(momentum_and_energy) = 0;
		std::get<1>(momentum_and_energy) = 0;
		std::get<2>(momentum_and_energy) = 0;
		for (auto& ball : physics.balls)
		{
			std::get<0>(momentum_and_energy) += ball.mass * ball.vx;
			std::get<1>(momentum_and_energy) += ball.mass * ball.vy;
			std::get<2>(momentum_and_energy) += 1.0f / 2.0f * ball.mass * (ball.vx * ball.vx + ball.vy * ball.vy);
		}
		return momentum_and_energy;
	}

	std::tuple<int, int> W_to_S(std::tuple<float, float> wCoord) {
		std::tuple<int, int> sCoord;
		std::get<0>(sCoord) = int((std::get<0>(wCoord) / (RIGHT_EDGE - LEFT_EDGE)) * ScreenWidth());
		std::get<1>(sCoord) = int((std::get<1>(wCoord) / (BOTTOM_EDGE - TOP_EDGE)) * ScreenHeight());
		return sCoord;
	}

	std::tuple<float, float> S_to_W(std::tuple<int, int> sCoord) {
		std::tuple<float, float> wCoord;
		std::get<0>(wCoord) = (std::get<0>(sCoord) * (RIGHT_EDGE - LEFT_EDGE)) / ScreenWidth();
		std::get<1>(wCoord) = (std::get<1>(sCoord) * (BOTTOM_EDGE - TOP_EDGE)) / ScreenHeight();
		return wCoord;
	}

public:
	bool OnUserCreate() override
	{
		srand(static_cast<unsigned int>(time(nullptr)));
		pixelsPerUnit = ScreenWidth() / (RIGHT_EDGE - LEFT_EDGE);
		sAppName = "Dropping Balls";
		spawnMainMenuObjects();
		return true;
	}

	bool OnUserUpdate(float fElapsedTime) override
	{
		Clear(olc::BLACK);

		//Handle the ball movements and collision physics
		if (!physics.timeStop)
		{
			auto start = std::chrono::system_clock::now();
			float fActualElapsedTime = physics.physicsHandling(fElapsedTime, true);
			auto end = std::chrono::system_clock::now();
			float elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<float>>(end - start).count();
			//Ratio signals how many seconds of sim time can be generated from 1 second of computation time
			float ratio = fActualElapsedTime / elapsed_seconds;
			std::this_thread::sleep_for(std::chrono::microseconds(20000 - (int)(fActualElapsedTime*1000000)));
			simulationTime += fActualElapsedTime;
		}
				
		//Cleanup routines
		physics.cleanup();

		//Handle the spawning of many balls
		massBallSpawn(fElapsedTime);

		//Input handling for main menu
		if (appState == AppState::STATE_MENU)
		{
			if ((GetMouseY() > ScreenHeight() / 2 - 200) && (GetMouseY() < ScreenHeight() / 2 - 190) && (GetMouse(0).bPressed))
			{
				despawnBalls();
				despawnObstacles();
				simulationTime = 0.0f;
				physics.score = 0;
				physics.gravity = 10.0f; //meters per second
				physics.restitutionFactor = 0.99f; //controls elasticity of collisions
				physics.drag = 0.8f;
				simulationTime = 0.0f;
				timeBetweenBalls = 0.05f;
				physics.ballVballCollisions = false;
				physics.invincibleObstacles = false;
				appState = AppState::STATE_GAME;
				gameState = GameState::GAMESTATE_NEWOBSTACLES;
				game_spawnNewObstacles();
				seconds_moved = 0.0f;
				spawnFrame();
			}

			if ((GetMouseY() > ScreenHeight() / 2 - 190) && (GetMouseY() < ScreenHeight() / 2 - 180) && (GetMouse(0).bPressed))
			{
				despawnBalls();
				despawnObstacles();
				timeBetweenBalls = 0.05f;
				physics.ballVballCollisions = true;
				physics.invincibleObstacles = true;
				physics.gravity = 10.0f; //meters per second
				physics.restitutionFactor = 0.99f; //controls elasticity of collisions
				physics.drag = 0.8f;
				simulationTime = 0.0f;
				appState = AppState::STATE_SANDBOX;
			}
		}
		//Input handling for sandbox mode
		else if (appState == AppState::STATE_SANDBOX)
		{
			//D: switch debug mode
			if (GetKey(olc::Key::D).bPressed)
			{
				debugMenu = !debugMenu;
			}

			//Debug Menu
			configurationInputHandling();

			spawnInputHandling();
		}
		//Input handling and logic for the actual game
		else if (appState == AppState::STATE_GAME)
		{
			if (gameState == GameState::GAMESTATE_END)
			{

			}
			else if (gameState == GameState::GAMESTATE_NEWOBSTACLES)
			{
				seconds_moved += fElapsedTime;
				if (seconds_moved <= moveUpInSeconds)
					game_moveObstacles(0.0f, Y_move * (fElapsedTime / moveUpInSeconds));
				else
				{
					game_moveObstacles(0.0f, Y_move * (((moveUpInSeconds - seconds_moved - fElapsedTime) / moveUpInSeconds)));
					gameState = GameState::GAMESTATE_AIM;
					for (auto& obstacle : physics.obstacles)
					{
						if (obstacle->type == ObstacleType::OBS_BALL)
							if (static_cast<BallObstacle*>(obstacle)->py > TOP_EDGE)
								gameState = GameState::GAMESTATE_END;

						if ((obstacle->type == ObstacleType::OBS_LINE) && !(static_cast<LineObstacle*>(obstacle)->isFrame))
							if ((static_cast<LineObstacle*>(obstacle)->sy > TOP_EDGE) || (static_cast<LineObstacle*>(obstacle)->ey > TOP_EDGE))
								gameState = GameState::GAMESTATE_END;
					}
				}
			}
			else if (gameState == GameState::GAMESTATE_AIM)
			{
				if (GetMouse(0).bPressed)
				{
					gameState = GameState::GAMESTATE_SHOOTING;
					ballsToSpawn = ballsPerShot;
					std::tuple<int, int> sCoord;
					std::get<0>(sCoord) = GetMouseX();
					std::get<1>(sCoord) = GetMouseY();
					targetCoord = S_to_W(sCoord);
					round++;
				}
				//Displays a prediction of the achieved score
				if (GetKey(olc::S).bHeld)
				{
					AI.physics.copyFrom(physics);
					std::tuple<float, float> wCoord;
					std::tuple<int, int> sCoord;
					std::get<0>(sCoord) = GetMouseX();
					std::get<1>(sCoord) = GetMouseY();
					wCoord = S_to_W(sCoord);
					AI.simulateShot(std::get<0>(wCoord), std::get<1>(wCoord), ballsPerShot, timeBetweenBalls);
					DrawString(ScreenWidth() / 2 - 30, ScreenHeight() / 2 - 200, std::to_string(AI.physics.score), olc::WHITE, 1);
				}
			}
			else if (gameState == GameState::GAMESTATE_SHOOTING)
			{
				if (physics.balls.size() == 0)
				{
					gameState = GameState::GAMESTATE_AIM;
					if ((round % spawnAfterRounds == 0) || (game_countObstacles() == 0))
					{
						gameState = GameState::GAMESTATE_NEWOBSTACLES;
						game_spawnNewObstacles();
						seconds_moved = 0.0f;
					}
				}
			}
			if (GetKey(olc::ESCAPE).bPressed)
			{
				returnToMenu();
			}
		}
		
		//DRAWING ROUTINE starts here

		//Draw balls and obstacles
		drawObjects();
		//Drawing routine for main menu
		if (appState == AppState::STATE_MENU)
		{
			DrawString(ScreenWidth() / 2 - 100, 40, "DROPPING BALLS", olc::WHITE, 2);
			DrawString(ScreenWidth() / 2 - 60, 60, "Author: Baran Oener", olc::WHITE, 1);
			if ((GetMouseY() > ScreenHeight() / 2 - 200) && (GetMouseY() < ScreenHeight() / 2 - 190))
				DrawString(ScreenWidth() / 2 - 30, ScreenHeight() / 2 - 200, "New Game", olc::GREEN, 1);
			else
				DrawString(ScreenWidth() / 2 - 30, ScreenHeight() / 2 - 200, "New Game", olc::WHITE, 1);
			if ((GetMouseY() > ScreenHeight() / 2 - 190) && (GetMouseY() < ScreenHeight() / 2 - 180))
				DrawString(ScreenWidth() / 2 - 50, ScreenHeight() / 2 - 190, "Sandbox Mode", olc::GREEN, 1);
			else
				DrawString(ScreenWidth() / 2 - 50, ScreenHeight() / 2 - 190, "Sandbox Mode", olc::WHITE, 1);
		}
		//Drawing routine for sandbox mode
		else if (appState == AppState::STATE_SANDBOX)
		{
			//SIMULATE the path of a shot ball
			std::tuple<float, float> wCoord;
			std::tuple<int, int> sCoord;
			std::get<0>(sCoord) = GetMouseX();
			std::get<1>(sCoord) = GetMouseY();
			wCoord = S_to_W(sCoord);	
			drawSimulatedBallPath(physics.spawn_x, physics.spawn_y, (std::get<0>(wCoord) - physics.spawn_x) * 2, (std::get<1>(wCoord) - physics.spawn_y) * 2, true, olc::WHITE);

			//Draw Debug Overlay
			if (debugMenu)
			{
				drawConfigurationMenu();
			}
		}

		//Draw routine for game
		else if (appState == AppState::STATE_GAME)
		{
			DrawString(ScreenWidth() - 100,  10, "Score: " + std::to_string(physics.score), olc::WHITE, 1);
			DrawString(10, 10, "Press ESC to Exit", olc::WHITE, 1);
			if (gameState == GameState::GAMESTATE_AIM)
			{
				//SIMULATE the path of a shot ball
				std::tuple<float, float> wCoord;
				std::tuple<int, int> sCoord;
				std::get<0>(sCoord) = GetMouseX();
				std::get<1>(sCoord) = GetMouseY();
				wCoord = S_to_W(sCoord);
				drawSimulatedBallPath(physics.spawn_x, physics.spawn_y, (std::get<0>(wCoord) - physics.spawn_x) * 2, (std::get<1>(wCoord) - physics.spawn_y) * 2, false, olc::WHITE);
			}
			else if (gameState == GameState::GAMESTATE_END)
			{
				DrawString(ScreenWidth() / 2 - 100, ScreenHeight() / 2, "GAME OVER", olc::WHITE, 2);
			}
		}

		return true;
	}

};

int main()
{
	BallsDrop game;
	if (game.Construct(600, 750, 1, 1))
	{
		game.Start();
	}
	return 0;
}