import numpy as np
from dataclasses import dataclass, field
import matplotlib.pyplot as plt

@dataclass
class Waypoint:
    location: np.ndarray
    velocity: np.ndarray = None
    acceleration: np.ndarray = None

    def checkIfDerivativesActive(self):
        return (self.checkIfVelocityActive() or self.checkIfAccelerationActive())
        
    def checkIfVelocityActive(self):
        return (self.velocity is not None)
    
    def checkIfAccelerationActive(self):
        return (self.acceleration is not None)
    
    def checkIfZeroVel(self):
        if self.velocity is not None and np.linalg.norm(self.velocity) <= 0:
            return True
        else:
            return False
        
    def print_waypoint(self, waypoint_num):
        print("Waypoint ", waypoint_num, ":")
        print(" location: ", self.location.flatten())
        print(" velocity: ", self.velocity.flatten())
        if self.acceleration is not None:
            print(" acceleration: ", self.acceleration.flatten())


@dataclass
class WaypointData:
    start_waypoint: Waypoint = field(default_factory=lambda: Waypoint(location=np.array([])))
    end_waypoint: Waypoint = field(default_factory=lambda: Waypoint(location=np.array([])))
    intermediate_locations: np.ndarray = None

    def get_waypoint_locations(self):
        point_sequence = self.start_waypoint.location
        if self.intermediate_locations is not None:
            point_sequence = np.concatenate((point_sequence, self.intermediate_locations),1)
        point_sequence = np.concatenate((point_sequence, self.end_waypoint.location),1)
        return point_sequence
    
    def get_distance(self):
        start_pos = self.start_waypoint.location
        end_pos = self.end_waypoint.location
        distance = np.linalg.norm(end_pos-start_pos)
        return distance

    def set_location_data(self, waypoint_sequence, set_vel = True):
        self.start_waypoint.location = waypoint_sequence[:,0][:,None]
        self.end_waypoint.location = waypoint_sequence[:,-1][:,None]
        if set_vel:
            self.start_waypoint.velocity = (waypoint_sequence[:,1] - waypoint_sequence[:,0])[:,None]
            self.end_waypoint.velocity = (waypoint_sequence[:,-1] - waypoint_sequence[:,-2])[:,None]
        if np.shape(waypoint_sequence)[1] > 2:
            self.intermediate_locations = waypoint_sequence[:,1:-1]

    def get_num_intermediate_waypoints(self):
        if self.intermediate_locations is not None:
            return np.shape(self.intermediate_locations)[1]
        else:
            return 0
        
def plot2D_waypoints(waypoint_data: WaypointData, ax, dot_size = 0.1, arrow_scale = 1):
    locations = waypoint_data.get_waypoint_locations()
    ax.scatter(locations[0,:],locations[1,:],color="b")
    if waypoint_data.start_waypoint.checkIfVelocityActive():
        start_pos = waypoint_data.start_waypoint.location
        if not waypoint_data.start_waypoint.checkIfZeroVel():
            start_dir = waypoint_data.start_waypoint.velocity
            ax.arrow(start_pos.item(0), start_pos.item(1), 
                    start_dir.item(0)*arrow_scale, start_dir.item(1)*arrow_scale, color = "b",
                    head_width = 0.5*arrow_scale, head_length= 0.5*arrow_scale)
            print("arrow_scale:", arrow_scale)
    if waypoint_data.end_waypoint.checkIfVelocityActive():
        end_pos = waypoint_data.end_waypoint.location
        if not waypoint_data.end_waypoint.checkIfZeroVel():
            end_dir = waypoint_data.end_waypoint.velocity
            ax.arrow(end_pos.item(0), end_pos.item(1), 
                    end_dir.item(0)*arrow_scale, end_dir.item(1)*arrow_scale, color = "b",
                    head_width = 0.5*arrow_scale, head_length= 0.5*arrow_scale)
    if waypoint_data.intermediate_locations is not None:
        ax.scatter(waypoint_data.intermediate_locations[0,:], 
                   waypoint_data.intermediate_locations[1,:],color="b", s=dot_size)

def plot3D_waypoints(waypoint_data: WaypointData, ax,  arrow_scale=2, arrow_length = None):
    locations = waypoint_data.get_waypoint_locations()
    ax.scatter(locations[0,:],locations[1,:],locations[2,:],color="b")
    if arrow_length is None:
        distance = waypoint_data.get_distance()
        arrow_length = distance/10*arrow_scale
    if waypoint_data.start_waypoint.checkIfVelocityActive():
        start_pos = waypoint_data.start_waypoint.location
        start_vel = waypoint_data.start_waypoint.velocity
        ax.quiver(start_pos.item(0), start_pos.item(1), start_pos.item(2), 
                  start_vel.item(0), start_vel.item(1), start_vel.item(2), 
                  length=arrow_length, normalize=True)
    if waypoint_data.end_waypoint.checkIfVelocityActive():
        end_pos = waypoint_data.end_waypoint.location
        end_vel = waypoint_data.end_waypoint.velocity
        ax.quiver(end_pos.item(0), end_pos.item(1), end_pos.item(2), 
                  end_vel.item(0), end_vel.item(1), end_vel.item(2), 
                  length=arrow_length, normalize=True)
    if waypoint_data.intermediate_locations is not None:
        ax.scatter(waypoint_data.intermediate_locations[0,:], 
                   waypoint_data.intermediate_locations[1,:],
                   waypoint_data.intermediate_locations[2,:],color="b")