#### Camera Transformation Maya Python Script
#### Author: Jaewon Song(songjaewon@kaist.ac.kr)
#
#### How-to-Install
# 1) Numpy Install
# 1-1) Unzip numpy.1.9.2.7z
# 1-2) Copy 'numpy' folder to Maya Python package folder
#      ex) C:\Program Files\Autodesk\Maya2017\Python\Lib\site-packages   
# 2) Run the script
#
##### How-to-Use
# 0) Positioning your camera with markers at the start frame
# 1) Run the entire script at once
# 2) Run the command:
#    camera_transform(camera_name, [marker_list], [start_frame, end_frame])
# ex)camera_transform('flwCAM1', ['Right_T', 'Top_T', 'Bottom_T', 'Left_T'], [1, 3489])

import maya.cmds as mc
import numpy as np
from numpy import *
from math import sqrt

def rigid_transform_3D(A, B):
    assert len(A) == len(B)

    N = A.shape[0]; # total points
    centroid_A = mean(A, axis=0)
    centroid_B = mean(B, axis=0)
    
    # centre the points
    AA = A - tile(centroid_A, (N, 1))
    BB = B - tile(centroid_B, (N, 1))

    # dot is matrix multiplication for array
    H = np.dot(transpose(AA),BB)
    U, S, Vt = linalg.svd(H)

    R = np.dot(Vt.T,U.T)

    # special reflection case
    if linalg.det(R) < 0:
       print "Reflection detected"
       Vt[2,:] *= -1
       R = np.dot(Vt.T,U.T)

    t = np.dot(-R,centroid_A.T) + centroid_B.T

    return R, t

def compute_camera_transform_with_points(point_list, init_pos_list, prev_frame, cur_frame):
    num_points = len(point_list)    
    prev_pos_list = init_pos_list
    cur_pos_list = []
    mc.currentTime(cur_frame)
    for point in point_list:
        cur_pos_list.append(mc.xform(point, q=True, t=True, ws=True))
    A = np.array(prev_pos_list)
    B = np.array(cur_pos_list)
    R, t = rigid_transform_3D(A,B)   
    m = [R[0][0], R[1][0], R[2][0], 0.0, R[0][1], R[1][1], R[2][1], 0.0, R[0][2], R[1][2], R[2][2], 0.0,  t[0], t[1], t[2], 1.0]    
    return m

def camera_transform(camera_name, marker_list, frame_range):
    if mc.objExists('camera_transform_locator'):
        mc.delete('camera_transform_locator')
    camera_loc = mc.spaceLocator(n='camera_transform_locator')
    mc.currentTime(frame_range[0])
    init_pos_list = []
    for marker in marker_list:
        init_pos_list.append(mc.xform(marker, q=True, t=True, ws=True))
    
    for frame in range(frame_range[0], frame_range[1]):
        m = compute_camera_transform_with_points(marker_list, init_pos_list, frame_range[0], frame+1)  
        mc.xform(camera_loc, m=m, ws=True)
        mc.setKeyframe(camera_loc)
    mc.currentTime(frame_range[0])
    mc.parentConstraint(camera_loc, camera_name, mo=True)
    
    
#camera_transform('s2_1_c001', ['Right_T', 'Top_T', 'Bottom_T', 'Left_T'], [1, 5700])
