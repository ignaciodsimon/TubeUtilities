#!/bin/sh

gcc `pkg-config --cflags gtk+-3.0` -o graphicInterface graphicInterface.c `pkg-config --libs gtk+-3.0`
#gcc `pkg-config --cflags gtk+-3.0` -o graphicInterface2 graphicInterface2.c `pkg-config --libs gtk+-3.0`
#gcc `pkg-config --cflags gtk+-3.0` -o openDialog openDialog.c `pkg-config --libs gtk+-3.0`
