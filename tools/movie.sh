ffmpeg -f image2 -r 40 -i movie/bender_%03d.png -c:v libx264 -pix_fmt yuv420p movie/bender.mp4
