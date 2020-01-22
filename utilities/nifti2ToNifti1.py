import nibabel as nib
import sys, os, shutil
import glob as glob

workdir = sys.argv[1]



if os.path.isdir(workdir):
    os.chdir(workdir)
else:
    sys.exit(workdir + " is not a directory")

A = glob.glob('*.nii*')

print('Converting all nifti-s in '+workdir +' to Nifti-1')
print('\nRead '+str(len(A))+ ' files, starting conversion: \n')


for a in A:
    print(a)
    shutil.copy(a,'bu_'+a)
    im2 = nib.load(a)
    im1 = nib.Nifti1Image(im2.dataobj, im2.affine, im2.header)
    im1.to_filename(a)
    os.remove('bu_'+a)
print('\nFinished, converted '+str(len(A))+ ' files.')
