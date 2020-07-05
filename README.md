# Dropping-Balls
A simple game and sandbox environment based on simulating the collision responses of balls.

The program is based on OneLoneCoder's olc::PixelGameEngine as well as his videos on collision handling of balls with some deviations.

The goal is simple - the player shoots balls, and tries to maximize the score by hitting (and destroying) as many obstacles as possible each turn. After a turn, new obstacles spawn from below and the entire grid moves upwards; the game is over when any obstacle reaches the top (taking inspiration from Tetris).

Some steps to potentially be worked on:
- Implementing a more precise way of detecting collisions. In order to make the collision responses consistent (meaning that exactly coinciding configuration of balls, velocity vectors and obstacles results in the exact same response), the physics algorithm currently uses a fixed "resolution/accuracy" by fixing the time steps in which the physics system is updated (decoupled from the framerate). This is a suboptimal solution, for several reasons that may be expanded on in time.
- Adding an "AI" player.
- Adding a two player mode.
- And, of course, more efficient code.

![Menu](https://github.com/BaranCanOener/Dropping-Balls/blob/master/Menu.PNG)
![Sandbox](https://github.com/BaranCanOener/Dropping-Balls/blob/master/Sandbox%20Mode.png)
![Game](https://github.com/BaranCanOener/Dropping-Balls/blob/master/Game.png)
