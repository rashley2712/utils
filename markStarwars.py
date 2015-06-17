#!/usr/bin/python
import pygame
from sys import exit
import time

size = width, height = 250, 150

pygame.init()
pygame.display.init()
pygame.display.set_caption( "Episode 7" )

screen = pygame.display.set_mode(size)

myfont = pygame.font.SysFont( "sans", 36 )
myfont2 = pygame.font.SysFont( "sans", 24 )

SW = 1450396800

while True:
  for event in pygame.event.get():
    if event.type == pygame.QUIT:
      pygame.quit()
      exit(0)
  screen.fill( (0,0,0) )
  t_now = time.time()
  t_str = "%08d s" %( SW - t_now )
  t_str_y = "%010.6f d" %( float(SW-t_now)/86400. )
  label = myfont.render(t_str, 1, (255,255,255) )
  label2 = myfont2.render( "Time till SW: ep vii", 1, (255,255,0) )
  label3 = myfont.render(t_str_y, 1, (255,255,255) )
  tw,th = pygame.Surface.get_size( label )
  tw2,th2 = pygame.Surface.get_size( label2 )
  tw3,th3 = pygame.Surface.get_size( label3 )
  screen.blit( label,  ( (width-tw)/2, 45 ) )
  screen.blit( label2, ( (width-tw2)/2, 7 ) )
  screen.blit( label3, ( (width-tw3)/2, 90 ) )
  pygame.display.flip()
  pygame.time.wait(1000)


