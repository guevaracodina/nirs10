spm_orthviews('AddBlobs',handle,XYZ,Z,mat)
% Adds blobs from a pointlist to the image specified by the handle(s).
% handle   - image number to add blobs to
% XYZ      - blob voxel locations
% Z        - blob voxel intensities
% mat      - matrix from voxels to millimeters of blob.
% name     - a name for this blob
% This method only adds one set of blobs, and displays them using a
% split colour table.

%%%%% EN VOXELS %%%%%%%%%
xyz = V.mat\[NIRS.SrcPos_MNI(1,:)';1]

 spm_orthviews('AddBlobs',gca,NIRS.SrcPos_MNI(1,:),1,V.mat);
 
 spm_orthviews('AddColouredBlobs',gca,NIRS.SrcPos_MNI(1,:),1,V.mat,[1 0 0],'BLOOOOOOOOOOOB');
