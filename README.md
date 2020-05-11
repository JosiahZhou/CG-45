# CG-45

Matthias Kortleven          -   A1
Celal  Karakoc              -   B
Jesse Jansen                -   B
Suleiman Kulane             -   A2
Shije Zhou                  -   A1
Daniel van Den Akker        -   B

Task distributions and work done

    Week 1:
    Accelaration structure                  -   Jesse/Celal         -   WIP
    Progress bar                            -   Jesse               -   Done
    Bounding boxes                          -   Celal               -   Done
    Tree structure                          -   Celal               -   Done
    Reflection/refraction                   -   Matthias            -   Working but buggy
    Recursive ray tracing                   -   Matthias            -   Done
    Soft shadows                            -   Suleiman      -   WIP
    Hard shadows                            -   Suleiman      -   WIP
    Diffuse stuff                           -   Josiah              -   WIP
    Multithreading                          -   Daniel              -   WIP

Bug report:
    Reflection refraction rays sometimes don't intersect with the right triangle, this is a bug in the tree of bounding boxes where currently only the first box that finds an intersection is checked, even though a triangle in another box can be closer to the origin of the ray.
    Reflection and refraction rays are not traced into open space (this is desired behavior, but makes it hard to see results of refraction/reflection)

    Week 2:
    Fix bounding box tree problem           -   Celal
    Make simple scene/materials             -   Matthias
