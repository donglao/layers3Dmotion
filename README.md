# Extending Layered Models to 3D Motion

This is the demo code for 'Extending Layered Models to 3D Motion (ECCV 2018)'. It performs the optimization described in the paper and outputs pixel level motion segmentation results.

## Getting Started

Simply run demo.m in Ubuntu.

In this demo, the optical flow code from [1] is used to estimate motion between images and the edge detector from [2] is used to generate initial region. We choose them since they are easy to set up and can run in different environment (many thanks to the original authors!). Both may be replaced by any state-of-the-art methods.

## To Be Done

This demo fully reproduces the optimization method described in the paper (equation (3)-(6)). To perform motion segmentation in longer sequences, a faster/simplified version (using motion only) will be released.

## Citation

```
@inproceedings{lao2018extending,
  title={Extending Layered Models to 3D Motion},
  author={Lao, Dong and Sundaramoorthi, Ganesh},
  booktitle={Proceedings of the European Conference on Computer Vision (ECCV)},
  pages={435--451},
  year={2018}
}
```

## References

[1] Sun, D., Roth, S., Black, M.J.: Secrets of optical flow estimation and their principles. In: Computer Vision and Pattern Recognition (CVPR), 2010 IEEE Conference on, IEEE (2010) 2432-2439

[2] Dollar, P., Zitnick, C.L.: Fast edge detection using structured forests. IEEE transactions on pattern analysis and machine intelligence 37(8) (2015) 1558-1570
