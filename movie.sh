ffmpeg -f image2 -r 10 -i bender_%02d.png -c:v libx264 -pix_fmt yuv420p bender.mp4
