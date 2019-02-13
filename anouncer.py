import time
import pygame

def drum():

    pygame.mixer.init()
    pygame.mixer.music.load('/home/pieter/Documents/Scriptie/utilities/utilities/sounds/BaDumTss.mp3')
    pygame.mixer.music.play()
    
    time.sleep(2)

