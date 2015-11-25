/*
 *  This file is part of esynth.
 *
 *  esynth is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  esynth is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with esynth.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _ID_FACTORY_GUARD
#define _ID_FACTORY_GUARD

class IdFactory
{
  public:
    IdFactory();
    IdFactory(unsigned int dictatedMin);

    unsigned int getNextId() { return current++; }
    unsigned int min() const { return minId; } 
    void reset() { current = minId; }

  private:
    unsigned int minId;
    unsigned int current;
   
    const static unsigned int DEFAULT_ID = 0;
};

#endif
