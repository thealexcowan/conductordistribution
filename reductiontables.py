import sys
from sage.all import *
import cProfile
import math
import random


M_MAX = 10
z1 = ZZ(1)


class TypeData:
    def __init__(self, data):
        self.font_scale_factor = 1
        self.__dict__.update(data)
        self.name = data.get('name', None)
        self.ratio = self.discriminant / self.conductor

    def __hash__(self):
        if self.name:
            return self.name
        else:
            return tuple(self.__dict__.items())

    def dict_key(self):
        return (self.symbol, int(log(self.conductor, self.base())))

    def base(self):
        return self.p
    
    def font_size(self):
        if self.p == 2:
            return 1.3 * self.font_scale_factor * self.p**8/self.modulus
        elif self.p == 3:
            return 1.0 * self.font_scale_factor * self.p**6/self.modulus
        else:
            return NotImplementedError


###

def get_locations(p, name):
    if p == 2:
        locations = get_locations_2(name)
    elif p == 3:
        locations = get_locations_3(name)
    else:
        locations = get_locations_p(name)
    return locations


def get_locations_2(name):
    if name == '2^12 I1':
        locations = [(642, 426), (646, 430), (650, 422), (654, 418), (658, 446), (662, 442), (666, 434), (670, 438), (674, 410), (678, 414), (682, 406), (686, 402), (690, 386), (694, 390), (698, 394), (702, 398), (706, 462), (710, 458), (714, 450), (718, 454), (722, 478), (726, 474), (730, 466), (734, 470), (738, 486), (742, 482), (746, 494), (750, 490), (754, 502), (758, 498), (762, 510), (766, 506)]
    
    elif name == '2^12 I2':
        locations = [(641, 431), (643, 429), (645, 425), (647, 427), (649, 419), (651, 417), (653, 423), (655, 421), (657, 441), (659, 443), (661, 445), (663, 447), (665, 439), (667, 437), (669, 433), (671, 435), (673, 413), (675, 415), (677, 411), (679, 409), (681, 401), (683, 403), (685, 405), (687, 407), (689, 389), (691, 391), (693, 387), (695, 385), (697, 399), (699, 397), (701, 393), (703, 395), (705, 457), (707, 459), (709, 461), (711, 463), (713, 453), (715, 455), (717, 451), (719, 449), (721, 473), (723, 475), (725, 477), (727, 479), (729, 469), (731, 471), (733, 467), (735, 465), (737, 481), (739, 483), (741, 485), (743, 487), (745, 489), (747, 491), (749, 493), (751, 495), (753, 499), (755, 497), (757, 503), (759, 501), (761, 507), (763, 505), (765, 511), (767, 509)]
    
    elif name == '2^12 I3':
        locations = []
        #[(640, 429), (641, 428), (642, 431), (643, 430), (644, 426), (645, 427), (646, 425), (647, 424), (648, 417), (649, 416), (650, 419), (651, 418), (652, 420), (653, 421), (654, 422), (655, 423), (656, 442), (657, 443), (658, 441), (659, 440), (660, 447), (661, 446), (662, 444), (663, 445), (664, 436), (665, 437), (666, 438), (667, 439), (668, 435), (669, 434), (670, 432), (671, 433), (672, 415), (673, 414), (674, 412), (675, 413), (676, 409), (677, 408), (678, 411), (679, 410), (680, 402), (681, 403), (682, 401), (683, 400), (684, 407), (685, 406), (686, 404), (687, 405), (688, 391), (689, 390), (690, 388), (691, 389), (692, 385), (693, 384), (694, 387), (695, 386), (696, 396), (697, 397), (698, 398), (699, 399), (700, 395), (701, 394), (702, 392), (703, 393), (704, 458), (705, 459), (706, 457), (707, 456), (708, 462), (709, 463), (710, 461), (711, 460), (712, 455), (713, 454), (714, 452), (715, 453), (716, 448), (717, 449), (718, 450), (719, 451), (720, 475), (721, 474), (722, 472), (723, 473), (724, 479), (725, 478), (726, 476), (727, 477), (728, 470), (729, 471), (730, 469), (731, 468), (732, 465), (733, 464), (734, 467), (735, 466), (736, 483), (737, 482), (738, 480), (739, 481), (740, 487), (741, 486), (742, 484), (743, 485), (744, 491), (745, 490), (746, 488), (747, 489), (748, 495), (749, 494), (750, 492), (751, 493), (752, 497), (753, 496), (754, 499), (755, 498), (756, 501), (757, 500), (758, 503), (759, 502), (760, 505), (761, 504), (762, 507), (763, 506), (764, 509), (765, 508), (766, 511), (767, 510)] # This is all the points

    elif name == 'I4*(4)':
        locations = [(644, 348), (652, 340), (660, 332), (668, 324), (676, 364), (684, 356), (692, 372), (700, 380), (708, 316), (716, 308), (724, 300), (732, 292), (740, 276), (748, 284), (756, 260), (764, 268)]

    elif name == 'I4*(6)':
        locations = [(144,48)]

    elif name == 'I5*(4)':
        locations = [(642, 342), (646, 338), (650, 346), (654, 350), (658, 322), (662, 326), (666, 334), (670, 330), (674, 358), (678, 354), (682, 362), (686, 366), (690, 382), (694, 378), (698, 374), (702, 370), (706, 306), (710, 310), (714, 318), (718, 314), (722, 290), (726, 294), (730, 302), (734, 298), (738, 282), (742, 286), (746, 274), (750, 278), (754, 266), (758, 270), (762, 258), (766, 262)]

    elif name == 'I5*(6)':
        locations = [(184,48), (168,60), (168,36)]
    
    elif name == 'I6*(4)':
        locations = [(641, 337), (643, 339), (645, 343), (647, 341), (649, 349), (651, 351), (653, 345), (655, 347), (657, 327), (659, 325), (661, 323), (663, 321), (665, 329), (667, 331), (669, 335), (671, 333), (673, 355), (675, 353), (677, 357), (679, 359), (681, 367), (683, 365), (685, 363), (687, 361), (689, 379), (691, 377), (693, 381), (695, 383), (697, 369), (699, 371), (701, 375), (703, 373), (705, 311), (707, 309), (709, 307), (711, 305), (713, 315), (715, 313), (717, 317), (719, 319), (721, 295), (723, 293), (725, 291), (727, 289), (729, 299), (731, 297), (733, 301), (735, 303), (737, 287), (739, 285), (741, 283), (743, 281), (745, 279), (747, 277), (749, 275), (751, 273), (753, 269), (755, 271), (757, 265), (759, 267), (761, 261), (763, 263), (765, 257), (767, 259)] # top right corners of 2x2 squares

    elif name == 'I6*(6)':
        locations = [(164, 48), (172, 42), (172, 54), (180, 34), (180, 62), (188, 38), (188, 58)]

    elif name == 'I7*(4)':
        locations = []
        #[(640, 338), (641, 339), (642, 336), (643, 337), (644, 341), (645, 340), (646, 342), (647, 343), (648, 350), (649, 351), (650, 348), (651, 349), (652, 347), (653, 346), (654, 345), (655, 344), (656, 325), (657, 324), (658, 326), (659, 327), (660, 320), (661, 321), (662, 323), (663, 322), (664, 331), (665, 330), (666, 329), (667, 328), (668, 332), (669, 333), (670, 335), (671, 334), (672, 352), (673, 353), (674, 355), (675, 354), (676, 358), (677, 359), (678, 356), (679, 357), (680, 365), (681, 364), (682, 366), (683, 367), (684, 360), (685, 361), (686, 363), (687, 362), (688, 376), (689, 377), (690, 379), (691, 378), (692, 382), (693, 383), (694, 380), (695, 381), (696, 371), (697, 370), (698, 369), (699, 368), (700, 372), (701, 373), (702, 375), (703, 374), (704, 309), (705, 308), (706, 310), (707, 311), (708, 305), (709, 304), (710, 306), (711, 307), (712, 312), (713, 313), (714, 315), (715, 314), (716, 319), (717, 318), (718, 317), (719, 316), (720, 292), (721, 293), (722, 295), (723, 294), (724, 288), (725, 289), (726, 291), (727, 290), (728, 297), (729, 296), (730, 298), (731, 299), (732, 302), (733, 303), (734, 300), (735, 301), (736, 284), (737, 285), (738, 287), (739, 286), (740, 280), (741, 281), (742, 283), (743, 282), (744, 276), (745, 277), (746, 279), (747, 278), (748, 272), (749, 273), (750, 275), (751, 274), (752, 270), (753, 271), (754, 268), (755, 269), (756, 266), (757, 267), (758, 264), (759, 265), (760, 262), (761, 263), (762, 260), (763, 261), (764, 258), (765, 259), (766, 256), (767, 257)] # all of them

    elif name == 'I7*(6)':
        locations = [(162, 41), (162, 55), (166, 43), (166, 53), (170, 48), (174, 45), (174, 51), (178, 37), (178, 59), (182, 39), (182, 57), (186, 33), (186, 63), (190, 35), (190, 61)]

    elif name == 'I8*(6)':
        locations = [(175 - z1/2, 48 - z1/2)] # because this is a 2x2 square
    
    else:
        locations = []
    
    return locations

    
def get_locations_3(name):
    if name == 'I1*':
        locations = [((2119+2146)/ZZ(2), (4360+4387)/ZZ(2)), (2038, 4468), (2092, 4522), (1957, 4549), (2011, 4603), (1822, 4657), (1741, 4738), (1633, 4873), (1687, 4927), (1579, 4981), (1498, 5062), (1903, 4819), (2038, 4279), (2092, 4225), (1957, 4198), (2011, 4144), (1822, 4090), (1741, 4009), (1633, 3874), (1687, 3820), (1579, 3766), (1498, 3685), (1903, 3928)]

    elif name == 'I2*':
        locations = [((2164+2173)/ZZ(2), (4369+4378)/ZZ(2)), (2065, 4252), (1984, 4171), (1840, 4072), (1759, 4018), (1876, 3928), (1660, 3847), (1597, 3775), (1507, 3703), (2065,4495), (1984, 4576), (1840, 4675), (1759, 4729), #(1876, 4918),
                     (1660, 4900), (1597, 4972), (1507, 5044), (1903,4792)]

    elif name == 'I3*':
        locations = [(2133 + 36 + 12 - z1/2, 4374 - z1/2)]

    elif name == '3^12 I1':
        locations = [(236, 161 + z1/2)]#[(235, 161 + z1/2)] # (175, 139), (202, 151), (202, 172), (175, 184)],

    else:
        locations = []

    return locations


def get_locations_p(name):
    locations = []
    return locations


###

def get_type_data(m_max=M_MAX, with_divisible=True):
    type_list = []
    type_list += get_type_data_2(m_max=with_divisible, with_divisible=with_divisible)
    type_list += get_type_data_3(m_max=with_divisible, with_divisible=with_divisible)
    type_list += get_type_data_p(m_max=with_divisible, with_divisible=with_divisible)
    return type_list


def get_type_data_2(m_max=M_MAX, with_divisible=True):
    type_list = []
    p = ZZ(2)
    
    # II
    td = {'name': 'II(4)',
          'symbol': 'II',
          'latex_name': r'$\mathrm{II}(4)$',
          'conductor': p**4,
          'discriminant': p**4,
          'modulus': p**2,
          'density': p**(-2), # v(b) = 0, v(disc) = 4, b - b^2 + a^2 - a^3 != 0 mod 4
          'locations': [(640,640), (256,896), (896,896)],
          }
    type_list.append(TypeData(td))
    td = {'name': 'II(6)',
          'symbol': 'II',
          'latex_name': r'$\mathrm{II}(6)$',
          'conductor': p**6,
          'discriminant': p**6,
          'modulus': p**2,
          'density': 3*p**(-4), # v(a) = 0, v(b) >= 1, b - b^2 + a^2 - a^3 != 0 mod 4  or  v(a) >= 1, v(b) = 1
          'locations': [(256,384), (640,128)],
          }
    type_list.append(TypeData(td))
    td = {'name': 'II(7)',
          'symbol': 'II',
          'latex_name': r'$\mathrm{II}(7)$',
          'conductor': p**7,
          'discriminant': p**7,
          'modulus': p**2,
          'density': p**(-4), # v(a) = 0, v(b) = 1, v(disc) = 7
          'locations': [(896,384)],
          }
    type_list.append(TypeData(td))
        
    # III
    td = {'name': 'III(3)',
          'symbol': 'III',
          'latex_name': r'$\mathrm{III}(3)$',
          'conductor': p**3,
          'discriminant': p**4,
          'modulus': p**2,
          'density': p**(-3), # v(a) = 0, v(b) = 0, b - b^2 + a^2 - a^3 == 0 mod 4, psi_3(a) != 0 mod 8  (# psi_3(x) = 3x^4 + 6ax^2 + 12bx - a^2, the 3-division polynomial)
                              # or  v(a) = 1, v(b) = 0, psi_3(a) != 0 mod 8
          'locations': [(384,640), (640,896)],
          }
    type_list.append(TypeData(td))
    td = {'name': 'III(5)',
          'symbol': 'III',
          'latex_name': r'$\mathrm{III}(5)$',
          'conductor': p**5,
          'discriminant': p**6,
          'modulus': p**2,
          'density': p**(-4), # v(a) = 0, v(b) >= 1, v(disc) = 6, b - b^2 + a^2 - a^3 == 0 mod 4
          'locations': [(896,128)],
          }
    type_list.append(TypeData(td))
    td = {'name': 'III(7)',
          'symbol': 'III',
          'latex_name': r'$\mathrm{III}(7)$',
          'conductor': p**7,
          'discriminant': p**8,
          'modulus': p**3,
          'density': p**(-5), # v(a) = 1, v(b) = 2
          'locations': [(384,192)],
          }
    type_list.append(TypeData(td))
    td = {'name': 'III(8)',
          'symbol': 'III',
          'latex_name': r'$\mathrm{III}(8)$',
          'conductor': p**8,
          'discriminant': p**9,
          'modulus': p**3,
          'density': p**(-5), # v(a) = 1, v(b) >= 3
          'locations': [(384,64)],
          }
    type_list.append(TypeData(td))

    # IV
    td = {'name': 'IV',
          'symbol': 'IV',
          'latex_name': r'$\mathrm{IV}$',
          'conductor': p**2,
          'discriminant': p**4,
          'modulus': p**2,
          'density': p**(-3), # v(a) = 0, v(b) = 0, b - b^2 + a^2 - a^3 != 0 mod 4, psi_3(a) != 0 mod 8
                              # or  v(a) >= 2, v(b) = 0,  b - b^2 + a^2 - a^3 != 0 mod 4
          'locations': [(128,640), (896,640)],
          }
    type_list.append(TypeData(td))

    # I0*
    td = {'name': 'I0*(4)',
          'symbol': 'I0*',
          'latex_name': r'$\mathrm{I}^{\!*}_0\!(4)$',
          'conductor': p**4,
          'discriminant': p**8,
          'modulus': p**4,
          'density': p**(-5), # v(a) = 0, v(b) = 1, v(disc) = 8, no r,t such that psi_3(r) = 0 mod 32 and psi_2(r,t) = 0 mod 64
                              # or  v(a) >= 2, v(b) = 2, no r,t such that psi_3(r) = 0 mod 32 and psi_2(r) - t^2 = 0 mod 64  (psi_2(x) = 4x^3 + 4ax + 4b, the 2-division polynomial)
          'locations': [(128,224),
                        (544,352), (544,480), (608,288), (608,416)],
          }
    type_list.append(TypeData(td))
    td = {'name': 'I0*(5)',
          'symbol': 'I0*',
          'latex_name': r'$\mathrm{I}^{\!*}_0\!(5)$',
          'conductor': p**5,
          'discriminant': p**9,
          'modulus': p**4,
          'density': p**(-6), # v(a) = 0, v(b) = 1, v(disc) = 9
          'locations': [(672,288), (736,384), (672,480)],
          }
    type_list.append(TypeData(td))
    td = {'name': 'I0*(6)',
          'symbol': 'I0*',
          'latex_name': r'$\mathrm{I}^{\!*}_0\!(6)$',
          'conductor': p**6,
          'discriminant': p**10,
          'modulus': p**4,
          'density': p**(-6), # v(a) >= 2, v(b) = 3
          'locations': [(128,96)],
          }
    type_list.append(TypeData(td))
    
    # I1*
    td = {'name': 'I1*',
          'symbol': 'I1*',
          'latex_name': r'$\mathrm{I}^{\!*}_1$',
          'conductor': p**3,
          'discriminant': p**8,
          'modulus': p**4,
          'density': p**(-6), # v(a) = 0, v(b) = 1, v(disc) = 8, there exists r,t such that psi_3(r) = 0 mod 32 and psi_2(r,t) = 0 mod 64, there exists r = 1 or 2 mod 4 s.t. psi_3(r) = 0 mod 32
                              # or  v(a) = 0, v(b) = 2, v(disc) = 8, there exists r = 1 or 2 mod 4 s.t. psi_3(r) = 0 mod 32
          'locations': [(192,160), (544,288), (608,352)],
          }
    type_list.append(TypeData(td))

    # I2*
    td = {'name': 'I2*(4)',
          'symbol': 'I2*',
          'latex_name': r'$\mathrm{I}^{\!*}_2\!(4)$',
          'conductor': p**4,
          'discriminant': p**10,
          'modulus': p**5,
          'density': p**(-8), # v(a) = 0, v(b) = 1, v(disc) = 6, there exists r = 1 or 2 mod 4 s.t. psi_3(r) = 0 mod 32
          'locations': [(752,304), (720,272), (688,336), (656,368)],
          }
    type_list.append(TypeData(td))
    td = {'name': 'I2*(6)',
          'symbol': 'I2*',
          'latex_name': r'$\mathrm{I}^{\!*}_2\!(6)$',
          'conductor': p**6,
          'discriminant': p**12,
          'modulus': p**5,
          'density': p**(-9), # v(a) = 2, v(b) >= 4, v(disc) = 12, there exists r = 1 or 2 mod 4 s.t. psi_3(r) = 0 mod 32  (v(b) >= 5 in fact)
          'locations': [(224,16)],
          }
    type_list.append(TypeData(td))
    td = {'name': 'I2*(7)',
          'symbol': 'I2*',
          'latex_name': r'$\mathrm{I}^{\!*}_2\!(7)$',
          'conductor': p**7,
          'discriminant': p**13,
          'modulus': p**5,
          'density': p**(-9), # v(a) = 2, v(b) = 4, v(disc) = 13
          'locations': [(224,48)],
          }
    type_list.append(TypeData(td))

    # I3*
    td = {'name': 'I3*(4)',
          'symbol': 'I3*',
          'latex_name': r'$\mathrm{I}^{\!*}_3\!(4)$',
          'conductor': p**4,
          'discriminant': p**11,
          'modulus': p**6,
          'density': p**(-9), # v(a) = 0, v(b) = 1, v(disc) = 11, there exists r = 1 or 2 mod 4 s.t. psi_3(r) = 0 mod 32
          'locations': [(760,280), (744,264), (728,312), (712,296),
                        (648,328), (664,344), (680,376), (696,360)],
          }
    type_list.append(TypeData(td))
    td = {'name': 'I3*(5)',
          'symbol': 'I3*',
          'latex_name': r'$\mathrm{I}^{\!*}_3\!(5)$',
          'conductor': p**5,
          'discriminant': p**12,
          'modulus': p**5,
          'density': p**(-9), # v(a) = 2, v(b) >= 4, v(disc) = 12, there exists r = 1 or 2 mod 4 s.t. psi_3(r) = 0 mod 32
          'locations': [(160,16)],
          }
    type_list.append(TypeData(td))

    # I4*
    td = {'name': 'I4*(4)',
          'symbol': 'I4*',
          'latex_name': r'$\mathrm{I}^{\!*}_4\!(4)$',
          'conductor': p**4,
          'discriminant': p**12,
          'modulus': p**7,
          'density': p**(-10), # v(a) = 0, v(b) = 1, v(disc) = 8+m, there exists r = 1 or 2 mod 4 s.t. psi_3(r) = 0 mod 32
          'locations': get_locations(2, 'I4*(4)'),
          }
    type_list.append(TypeData(td))
    td = {'name': 'I4*(6)',
          'symbol': 'I4*',
          'latex_name': r'$\mathrm{I}^{\!*}_4\!(6)$',
          'conductor': p**6,
          'discriminant': p**14,
          'modulus': p**5,
          'density': p**(-10), # v(a) = 2, v(b) = 4, v(disc) = 10+m
          'locations': get_locations(2, 'I4*(6)'),
          }
    type_list.append(TypeData(td))
    
    # Im*
    for m in range(5, m_max+1):
        td = {'name': 'I' + str(m) + '*(4)',
              'symbol': 'I' + str(m) + '*',
              'latex_name': r'$\mathrm{I}^{\!*}_' + str(m) + r'\!(4)$',
              'conductor': p**4,
              'discriminant': p**(8+m),
              'modulus': p**(3+m),
              'density': p**(-m-6), # v(a) = 0, v(b) = 1, v(disc) = 8+m, there exists r = 1 or 2 mod 4 s.t. psi_3(r) = 0 mod 32
              'locations': get_locations(2, 'I' + str(m) + '*(4)'),
              }
        type_list.append(TypeData(td))
        td = {'name': 'I' + str(m) + '*(6)',
              'symbol': 'I' + str(m) + '*',
              'latex_name': r'$\mathrm{I}^{\!*}_' + str(m) + r'\!(6)$',
              'conductor': p**6,
              'discriminant': p**(10+m),
              'modulus': p**(2+m),
              'density': p**(-m-6), # v(a) = 2, v(b) = 4, v(disc) = 10+m
              'locations': get_locations(2, 'I' + str(m) + '*(6)'),
              }
        type_list.append(TypeData(td))
    
    # IV*
    td = {'name': 'IV*',
          'symbol': 'IV*',
          'latex_name': r'$\mathrm{IV}^{\!*}$',
          'conductor': p**2,
          'discriminant': p**8,
          'modulus': p**4,
          'density': p**(-6), # TODO
          'locations': [(64,160), (544,416), (608,480)],
          }
    type_list.append(TypeData(td))

    # III*
    td = {'name': 'III*(3)',
          'symbol': 'III*',
          'latex_name': r'$\mathrm{III}^{\!*}\!(3)$',
          'conductor': p**3,
          'discriminant': p**10,
          'modulus': p**5,
          'density': p**(-8), # TODO
          'locations': [(656, 400), (688, 432), (720, 496), (752, 464)],
          }
    type_list.append(TypeData(td))
    td = {'name': 'III*(5)',
          'symbol': 'III*',
          'latex_name': r'$\mathrm{III}^{\!*}\!(5)$',
          'conductor': p**5,
          'discriminant': p**12,
          'modulus': p**5,
          'density': p**(-9), # TODO
          'locations': [(96,48)],
          }
    type_list.append(TypeData(td))
    td = {'name': 'III*(7)',
          'symbol': 'III*',
          'latex_name': r'$\mathrm{III}^{\!*}\!(7)$',
          'conductor': p**7,
          'discriminant': p**14,
          'modulus': p**6,
          'density': p**(-10), # TODO
          'locations': [(96,24)],
          }
    type_list.append(TypeData(td))
    td = {'name': 'III*(8)',
          'symbol': 'III*',
          'latex_name': r'$\mathrm{III}^{\!*}\!(8)$',
          'conductor': p**8,
          'discriminant': p**15,
          'modulus': p**6,
          'density': p**(-10), # TODO
          'locations': [(96,8)],
          }
    type_list.append(TypeData(td))

    # II*
    td = {'name': 'II*(3)',
          'symbol': 'II*',
          'latex_name': r'$\mathrm{II}^{\!*}\!(3)$',
          'conductor': p**3,
          'discriminant': p**11,
          'modulus': p**6,
          'density': p**(-9), # TODO
          'locations': [(648, 440), (664, 424), (680, 392), (696, 408), (712, 472), (728, 456), (744, 504), (760, 488)],
          }
    type_list.append(TypeData(td))
    td = {'name': 'II*(4)',
          'symbol': 'II*',
          'latex_name': r'$\mathrm{II}^{\!*}\!(4)$',
          'conductor': p**4,
          'discriminant': p**12,
          'modulus': p**6,
          'density': p**(-10), # TODO
          'locations': [(32,56)],
          }
    type_list.append(TypeData(td))
    td = {'name': 'II*(6)',
          'symbol': 'II*',
          'latex_name': r'$\mathrm{II}^{\!*}\!(6)$',
          'conductor': p**6,
          'discriminant': p**14,
          'modulus': p**6,
          'density': p**(-10), # TODO
          'locations': [(32,24)],
          }
    type_list.append(TypeData(td))
    
    # I0
    td = {'name': '2^12 I0',
          'symbol': 'I0',
          'latex_name': r'$2^{12}\mathrm{I}_0$',
          'conductor': z1,
          'discriminant': p**12,
          'modulus': p**7,
          'density': p**(-9), # v(a) = 0, v(b) = 1, v(disc) = 12  or  v(a) >= 4, v(b) = 4, v(disc) = 12  (modulus 2^6)
          'locations': [(32,40)] +
                       [(644,420), (652,428), (660,436), (668,444),
                        (676,404), (684,412), (692,396), (700,388),
                        (708,452), (716,460), (724,468), (732,476),
                        (740,492), (748,484), (756,508), (764,500)],
          'color': 'black',
          }
    type_list.append(TypeData(td))
    
    # Im
    for m in range(1, m_max+1):
        td = {'name': '2^12 I' + str(m),
              'symbol': 'I' + str(m),
              'latex_name': r'$2^{12}\mathrm{I}_' + str(m) + '$',
              'conductor': p,
              'discriminant': p**(12+m),
              'modulus': p**(7+m),
              'density': p**(-m-10), # v(a) = 0, v(b) = 1, v(disc) = 12+m
              'locations': get_locations(p, '2^12 I' + str(m)),
              }
        type_list.append(TypeData(td))

    if with_divisible:
        # Divisible
        td = {'name': 'divisible',
              'symbol': 'divisible',
              'latex_name': r'$2^4 \mid a \text{ and } 2^6 \mid b$',
              'conductor': -1, 
              'discriminant': -1,
              'modulus': p**6,
              'density': p**(-10), # v(a) >= 4, v(b) >= 6
              'locations': [(32,8)],
              'color': 'lightgrey',
              }
        type_list.append(TypeData(td))

    for td in type_list:
        td.p = p

    return type_list



def get_type_data_3(m_max=M_MAX, with_divisible=True):
    type_list = []
    p = ZZ(3)
    
    # I0
    td = {'name': 'I0',
          'symbol': 'I0',
          'latex_name': r'$\mathrm{I}_0$',
          'conductor': z1,
          'discriminant': z1,
          'modulus': p,
          'density': 2*p**(-1),
          'locations': [(2*3**7, (3**8-1)/2)],
          'color': 'black',
          }
    type_list.append(TypeData(td))
    # See also 3^12 I0 below
    
    # II
    td = {'name': 'II(3)',
          'symbol': 'II',
          'latex_name': r'$\mathrm{II}(3)$',
          'conductor': p**3,
          'discriminant': p**3,
          'modulus': p**2,
          'density': 4*p**(-3), # TODO
          'locations': [(2*3**6, 2*3**6), (3**6 - z1/2, 2*3**7 - z1/2), ((3**7-1)/2, 3**8 - (3**6+1)/2)],
          }
    type_list.append(TypeData(td))
    td = {'name': 'II(4)',
          'symbol': 'II',
          'latex_name': r'$\mathrm{II}(4)$',
          'conductor': p**4,
          'discriminant': p**4,
          'modulus': p**2,
          'density': 4*p**(-4),
          'locations': [(3**7 - (3**6+1)/2, 3**7 + 3**6), (3**7 - (3**6+1)/2, 3**8 - 3**6)],
          }
    type_list.append(TypeData(td))
    td = {'name': 'II(5)',
          'symbol': 'II',
          'latex_name': r'$\mathrm{II}(5)$',
          'conductor': p**5,
          'discriminant': p**5,
          'modulus': p**2,
          'density': 2*p**(-4),
          'locations': [((3**6-1)/2, 2*3**6 - z1/2)],
          }
    type_list.append(TypeData(td))

    # III
    td = {'name': 'III',
          'symbol': 'III',
          'latex_name': r'$\mathrm{III}$',
          'conductor': p**2,
          'discriminant': p**3,
          'modulus': p**2,
          'density': 2*p**(-3),
          'locations': [(2*3**6 - z1/2, (3**6-1)/2),
                        ((3**6-1)/2, 3**7 + (3**6-1)/2), (3**6 + (3**6-1)/2, 3**7 + 3**6 + (3**6-1)/2),
                        ((3**6-1)/2, 3**8 - (3**6+1)/2), (3**6 + (3**6-1)/2, 3**8 - 3**6 - (3**6+1)/2)],
          }
    type_list.append(TypeData(td))

    # IV
    td = {'name': 'IV(3)',
          'symbol': 'IV',
          'latex_name': r'$\mathrm{IV}(3)$',
          'conductor': p**3,
          'discriminant': p**5,
          'modulus': p**3,
          'density': 4*p**(-5),
          'locations': [(2*3**6 + 3**5 - z1/2, 2*3**7 - z1/2),
                        (3**7 - 3**5 - z1/2, 3**7 + 2*3**6 + (3**5+1)/2),
                        (3**7 - 3**5 - z1/2, 3**8 - 2*3**6 - (3**5+1)/2)],
          }
    type_list.append(TypeData(td))
    td = {'name': 'IV(4)',
          'symbol': 'IV',
          'latex_name': r'$\mathrm{IV}(4)$',
          'conductor': p**4,
          'discriminant': p**6,
          'modulus': p**3,
          'density': 4*p**(-6),
          'locations': [(2*3**5, 2*3**5)],
          }
    type_list.append(TypeData(td))
    td = {'name': 'IV(5)',
          'symbol': 'IV',
          'latex_name': r'$\mathrm{IV}(5)$',
          'conductor': p**5,
          'discriminant': p**7,
          'modulus': p**3,
          'density': 2*p**(-6),
          'locations': [((3**5-1)/2, 2*3**5)],
          }
    type_list.append(TypeData(td))

    # I0*
    td = {'name': 'I0*',
          'symbol': 'I0*',
          'latex_name': r'$\mathrm{I}^{\!*}_0$',
          'conductor': p**2,
          'discriminant': p**6,
          'modulus': p**4,
          'density': 2*p**(-5),
          # 'locations': [(2*3**5, (3**5-1)/2),
          #               (2*3**6 + 3**5 - (3**4+1)/2, 3**7 + 2*3**6 + (3**4-1)/2),
          #               (2*3**6 + (3**4-1)/2, 3**7 + 2*3**6 + 3**5 - (3**4+1)/2),
          #               (2*3**6 + 3**5 - (3**4+1)/2, 2*3**7 + 3**6 - (3**4+1)/2),
          #               (2*3**6 + (3**4-1)/2, 2*3**7 + 3**6 - 3**5 + (3**4-1)/2),
          #               (2*3**6 + (3**6-1)/2, 3**7 + 2*3**6 + (3**6-1)/2),
          #               (2*3**6 + (3**6-1)/2, 2*3**7 + 3**6 - (3**6+1)/2),
          #               (3**7 - 3**5 + 3**4 - z1/2, 2*3**7 - z1/2)],
          'locations': [(1498, 3847), (1660, 3685), (1822, 4009), ((1984+2065)/ZZ(2), (4333+4414)/ZZ(2)),
                        (1498, 4900), (1660, 5062), (1822, 4738), ((364+607)/ZZ(2), 121)],
          'font_scale_factor': 1.6,
          }
    type_list.append(TypeData(td))

    # Im*
    for m in range(1, m_max):
        td = {'name': 'I' + str(m) + '*',
              'symbol': 'I' + str(m) + '*',
              'latex_name': r'$\mathrm{I}^{\!*}_' + str(m) + '$',
              'conductor': p**2,
              'discriminant': p**(6+m),
              'modulus': p**(4+m),
              'density': 4*p**(-m-6),
              'locations': get_locations(3, 'I' + str(m) + '*'),
              'font_scale_factor': 1.6 if m in [1,2] else 1.0,
              }
        type_list.append(TypeData(td))
    
    # IV*
    td = {'name': 'IV*(3)',
          'symbol': 'IV*',
          'latex_name': r'$\mathrm{IV}^{\!*}\!(3)$',
          'conductor': p**3,
          'discriminant': p**9,
          'modulus': p**5,
          'density': 4*p**(-8),
          'locations': [(2*3**4, 2*3**3), (3**4, 2*3**4), (3**4 + (3**4-1)/2, 3**5 - (3**3+1)/2)],
          'font_scale_factor': 1.9,
          }
    type_list.append(TypeData(td))
    td = {'name': 'IV*(4)',
          'symbol': 'IV*',
          'latex_name': r'$\mathrm{IV}^{\!*}\!(4)$',
          'conductor': p**4,
          'discriminant': p**10,
          'modulus': p**5,
          'density': 4*p**(-9),
          'locations': [(3**5 - (3**4+1)/2, 3**4 + 3**3), (3**5 - (3**4+1)/2, 3**5 - 3**3)],
          'font_scale_factor': 1.9,
          }
    type_list.append(TypeData(td))
    td = {'name': 'IV*(5)',
          'symbol': 'IV*',
          'latex_name': r'$\mathrm{IV}^{\!*}\!(5)$',
          'conductor': p**5,
          'discriminant': p**11,
          'modulus': p**5,
          'density': 2*p**(-9),
          'locations': [((3**4 - 1)/2, 3**4 - 3**3)],
          'font_scale_factor': 1.9,
          }
    type_list.append(TypeData(td))

    # III*
    td = {'name': 'III*',
          'symbol': 'III*',
          'latex_name': r'$\mathrm{III}^{\!*}$',
          'conductor': p**2,
          'discriminant': p**9,
          'modulus': p**5,
          'density': 2*p**(-8),
          'locations': [(2*3**4, (3**3-1)/2),
                        ((3**4-1)/2, 3**4 + (3**3-1)/2), ((3**4-1)/2, 3**5 - (3**3+1)/2),
                        ((3**4-1)/2 + 3**4, 3**4 + 3**3 + (3**3-1)/2), ((3**4-1)/2 + 3**4, 3**5 - 3**3 - (3**3+1)/2)],
          'font_scale_factor': 1.9,
          }
    type_list.append(TypeData(td))

    # II*
    td = {'name': 'II*(3)',
          'symbol': 'II*',
          'latex_name': r'$\mathrm{II}^{\!*}\!(3)$',
          'conductor': p**3,
          'discriminant': p**11,
          'modulus': p**6,
          'density': 4*p**(-10),
          'locations': [(3**5 - 2*3**3, 2*3**4 - z1/2)],# (3**5 - (3**3-1)/2, 2*3**4 - (3**3+1)/2), (3**5 - (3**3-1)/2, 2*3**4 + (3**3-1)/2)],
          'font_scale_factor': 2.5,
          }
    type_list.append(TypeData(td))
    td = {'name': 'II*(4)',
          'symbol': 'II*',
          'latex_name': r'$\mathrm{II}^{\!*}\!(4)$',
          'conductor': p**4,
          'discriminant': p**12,
          'modulus': p**6,
          'density': 4*p**(-11),
          'locations': [(2*3**3, 2*3**2-1)],
          'font_scale_factor': 2.5,
          }
    type_list.append(TypeData(td))
    td = {'name': 'II*(5)',
          'symbol': 'II*',
          'latex_name': r'$\mathrm{II}^{\!*}\!(5)$',
          'conductor': p**5,
          'discriminant': p**13,
          'modulus': p**6,
          'density': 2*p**(-11),
          'locations': [(13,18-1)],
          'font_scale_factor': 2.5,
          }
    type_list.append(TypeData(td))

    # 3^12 I0
    td = {'name': '3^12 I0',
          'symbol': 'I0',
          'latex_name': r'$3^{12}\,\mathrm{I}_0$',
          'conductor': z1,
          'discriminant': p**12,
          'modulus': p**7,
          'density': 4*p**(-11),
          'locations': [(222, 162 - z1/2)],#[(220, 162 - z1/2)],
          'color': 'black',
          'font_scale_factor': 1.3,
          }
    type_list.append(TypeData(td))
    
    # Im
    for m in range(1, m_max+1):
        td = {'name': '3^12 I' + str(m),
              'symbol': 'I' + str(m),
              'latex_name': r'$3^{12}\,\mathrm{I}_' + str(m) + '$',
              'conductor': p,
              'discriminant': p**(12+m),
              'modulus': p**(7+m),
              'density': 4*p**(-m-11), # v(a) = 3, v(b) = 3, v(disc) = 12+m
              'locations': get_locations(3, '3^12 I' + str(m)),
              'font_scale_factor': 0.7,
              }
        type_list.append(TypeData(td))

    if with_divisible:
        # Divisible
        td = {'name': 'divisible',
              'symbol': 'divisible',
              'latex_name': r'$3^4 \mid a \text{ and } 3^6 \mid b$',
              'conductor': -1, 
              'discriminant': -1,
              'modulus': p**6,
              'density': p**(-10), # v(a) >= 4, v(b) >= 6
              'locations': [((3**4-1)/2, (3**2-1)/2)],
              'color': 'lightgrey',
              }
        type_list.append(TypeData(td))

    for td in type_list:
        td.p = p
    
    return type_list



def types_p(m_max=M_MAX, with_divisible=True, p=None):
    type_list = []
    if p is None:
        p = ZZ['p'].gen()
    else:
        p = ZZ(p)

    # I0
    td = {'name': 'I0',
          'symbol': 'I0',
          'latex_name': r'$\mathrm{I}_0$',
          'conductor': z1,
          'discriminant': z1,
          'modulus': p,
          'density': 1 - 1/p,
          'locations': [], # v(disc) = 0
          'color': 'black',
          }
    type_list.append(TypeData(td))

    # Im
    for m in range(1, m_max+1):
        td = {'name': 'I' + str(m),
              'symbol': 'I' + str(m),
              'latex_name': r'$\mathrm{I}_' + str(m) + '$',
              'conductor': p,
              'discriminant': p**m,
              'modulus': p**(1+m),
              'density': p**(-m)*(1 - 1/p)**2, # v(a) = v(b) = 0, v(disc) = m
              'locations': [], # TODO
              }
        type_list.append(TypeData(td))

    # II
    td = {'name': 'II',
          'symbol': 'II',
          'latex_name': r'$\mathrm{II}$',
          'conductor': p**2,
          'discriminant': p**2,
          'modulus': p**2,
          'density': p**(-2)*(1 - 1/p), # v(a) >= 1, v(b) = 1
          'locations': [], # TODO
          }
    type_list.append(TypeData(td))

    # III
    td = {'name': 'III',
          'symbol': 'III',
          'latex_name': r'$\mathrm{III}$',
          'conductor': p**2,
          'discriminant': p**3,
          'modulus': p**2,
          'density': p**(-3)*(1 - 1/p), # v(a) = 1, v(b) >= 2
          'locations': [], # TODO
          }
    type_list.append(TypeData(td))

    # IV
    td = {'name': 'IV',
          'symbol': 'IV',
          'latex_name': r'$\mathrm{IV}$',
          'conductor': p**2,
          'discriminant': p**4,
          'modulus': p**3,
          'density': p**(-4)*(1 - 1/p), # v(a) >= 2, v(b) = 2
          'locations': [], # TODO
          }
    type_list.append(TypeData(td))

    # I0*
    td = {'name': 'I0*',
          'symbol': 'I0*',
          'latex_name': r'$\mathrm{I}^{\!*}_0$',
          'conductor': p**2,
          'discriminant': p**6,
          'modulus': p**4,
          'density': p**(-5)*(1 - 1/p), # v(a) >= 2, v(b) >= 3, v(disc) = 6
          'locations': [], # TODO
          }
    type_list.append(TypeData(td))

    # Im*
    for m in range(1, m_max+1):
        td = {'name': 'I' + str(m) + '*',
              'symbol': 'I' + str(m) + '*',
              'latex_name': r'$\mathrm{I}^{\!*}_' + str(m) + '$',
              'conductor': p**2,
              'discriminant': p**(6+m),
              'modulus': p**(4+m),
              'density': p**(-m-5)*(1 - 1/p)**2, # v(a) = 2, v(b) = 3, v(disc) = 6+m
              'locations': [], # TODO
              }
        type_list.append(TypeData(td))

    # IV*
    td = {'name': 'IV*',
          'symbol': 'IV*',
          'latex_name': r'$\mathrm{IV}^{\!*}$',
          'conductor': p**2,
          'discriminant': p**8,
          'modulus': p**5,
          'density': p**(-7)*(1 - 1/p), # v(a) >= 3, v(b) = 4
          'locations': [], # TODO
          }
    type_list.append(TypeData(td))

    # III*
    td = {'name': 'III*',
          'symbol': 'III*',
          'latex_name': r'$\mathrm{III}^{\!*}$',
          'conductor': p**2,
          'discriminant': p**9,
          'modulus': p**5,
          'density': p**(-8)*(1 - 1/p), # v(a) = 3, v(b) >= 5
          'locations': [], # TODO
          }
    type_list.append(TypeData(td))

    # II*
    td = {'name': 'II*',
          'symbol': 'II*',
          'latex_name': r'$\mathrm{II}^{\!*}$',
          'conductor': p**2,
          'discriminant': p**10,
          'modulus': p**6,
          'density': p**(-9)*(1 - 1/p), # v(a) >= 4, v(b) = 5
          'locations': [], # TODO
          }
    type_list.append(TypeData(td))

    if with_divisible:
        # Divisible
        td = {'name': 'divisible',
              'symbol': 'divisible',
              'latex_name': r'$p^4\mid a \text{ and } p^6 \mid b$',
              'conductor': -1, 
              'discriminant': -1,
              'modulus': p**6,
              'density': p**(-10), # v(a) >= 4, v(b) >= 6
              'locations': [], # TODO
              'color': 'lightgrey',
              }
        type_list.append(TypeData(td))

    for td in type_list:
        td.p = p
    
    return type_list
