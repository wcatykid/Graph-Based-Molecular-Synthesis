/*
 *  This file is part of synth.
 *
 *  synth is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  synth is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with synth.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _TIMED_HASH_MAP_GUARD
#define _TIMED_HASH_MAP_GUARD 1


#include <utility>
#include <map>
#include <iostream>
#include <sstream>


#include "TimedLikeValueContainer.h"


class TimedHashMap
{
  private:
    TimedLikeValueContainer** like_containers;

    static unsigned SIZE;
    static unsigned rel_UPPERBOUND;
    static unsigned ABS_UPPERBOUND;
    unsigned int sz;
    unsigned int currentTime;

    //
    // Trusted hash function for strings (found online)
    //
    unsigned int hash(const char* str) const
    {
        unsigned long hash = 5381;
        int c;

        while (c = *str++)
        {
            hash = ((hash << 5) + hash) + c; /* hash * 33 + c */
        }

        return hash % SIZE;
    }

    //
    // Remove all the old values from the container
    //
    void PurgeAndRenumber()
    {
        std::cerr << "Purging" << std::endl;
    
        // Purge all the hashing containers
        for (int index = 0; index < SIZE; index++)
        {
            if (like_containers[index] != 0)
            {
                like_containers[index]->PurgeAndRenumber();

                // If this list is now empty, kill it.
                if (like_containers[index]->empty())
                {
                    delete like_containers[index];
                    like_containers[index] = 0;
                }
            }
        }

        // Update the size of the number of elements in this container to the relative max.
        sz = rel_UPPERBOUND;

        // Reset time to the sliding window upper-bound. We then subtract down to 0.
        currentTime = ABS_UPPERBOUND - rel_UPPERBOUND;
    }	

    //
    // Handle collisions with linear probing
    //
    bool privateAdd(const std::string& that)
    {
        //
        // Do we need to purge before addition?
        // (this updates sz)
        if (sz == ABS_UPPERBOUND) PurgeAndRenumber();
        
        //
        // Addition can be performed normally
        //
        unsigned int hashVal = hash(that.c_str()) % SIZE;

        // Do we need to create this list first?
        if (like_containers[hashVal] == 0)
        {
            like_containers[hashVal] = new TimedLikeValueContainer();
        }
        
        like_containers[hashVal]->add(that, currentTime);
        
        // Count down to 0.
        currentTime--;

        return true;
    }

	
  public:
    static void SetThresholds(unsigned relUB, unsigned absUB)
    {
        rel_UPPERBOUND = relUB;
        ABS_UPPERBOUND = absUB;
        SIZE = absUB;
        
        TimedLikeValueContainer::SetThresholds(relUB, absUB);
    }

    // Init to the absolute upper bound so we can populate the container completely.    
    TimedHashMap() : sz(0), currentTime(ABS_UPPERBOUND)
    {
        like_containers = new TimedLikeValueContainer*[SIZE];
		
        for (int i = 0; i < SIZE; i++)
        {
            like_containers[i] = 0;
        }
    }
	
    //
    // Delete all entries
    //
    ~TimedHashMap()
    {
        for (int i = 0; i < SIZE; i++)
        {
            if (like_containers[i]) delete like_containers[i];
        }

        delete[] like_containers;
    }

    unsigned int size() const { return sz; }

    bool add(const std::string& that)
    {
        bool added = privateAdd(that);

        if (added) sz++;

        return added;
    }

    //
    // Should only be called once
    //
    bool contains(const std::string& that) const
    {
        unsigned int hashVal = hash(that.c_str()) % SIZE;

        if (like_containers[hashVal] == 0) return false;
			
        // The list exists, inquire if the key matches our key.
        return like_containers[hashVal]->contains(that);
    }
    
    std::string toString() const
    {
        std::ostringstream oss;
        
        for (int index = 0; index < SIZE; index++)
        {
            if (like_containers[index] != 0)
            {
                oss << index << ": ";
                oss << *like_containers[index] << std::endl;
            }
        }
        
        return oss.str();
    }
    
    friend std::ostream& operator<<(std::ostream& os, const TimedHashMap& map)
    {
        os << map.toString();
        
        return os;
    }
};

#endif
