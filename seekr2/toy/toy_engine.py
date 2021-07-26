"""
toy_engine.py


"""

class BrownianParticle():
    def __init__(self, mass, dimensions, position, velocity, diffusion):
        self.mass = mass
        self.dimensions = dimensions
        assert len(position.shape) == dimensions
        self.position = position
        self.velocity = velocity
        self.diffusion = diffusion
        return