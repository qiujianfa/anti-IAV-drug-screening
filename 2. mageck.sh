#!/bin/bash

mageck test -k mageck_input.txt -t A_umicounts.p1,B_umicounts.p1,D_umicounts.p1 -c A_umicounts.p0 -n ./output/mageck_p1
mageck test -k mageck_input.txt -t A_umicounts.p2,B_umicounts.p2,D_umicounts.p2 -c A_umicounts.p0 -n ./output/mageck_p2
mageck test -k mageck_input.txt -t A_umicounts.p3,B_umicounts.p3,D_umicounts.p3 -c A_umicounts.p0 -n ./output/mageck_p3
mageck test -k mageck_input.txt -t A_umicounts.p4,B_umicounts.p4,D_umicounts.p4 -c A_umicounts.p0 -n ./output/mageck_p4
mageck test -k mageck_input.txt -t A_umicounts.p5,B_umicounts.p5,D_umicounts.p5 -c A_umicounts.p0 -n ./output/mageck_p5

