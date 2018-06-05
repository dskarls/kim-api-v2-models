   switch(GetComputeIndex(isComputeProcess_dEdr,
                          isComputeProcess_d2Edr2,
                          isComputeEnergy,
                          isComputeForces,
                          isComputeParticleEnergy))
   {
      case 0:
         ier = Compute< false, false,
                        false, false,
                        false >(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1:
         ier = Compute< false, false,
                        false, false,
                        true >(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 2:
         ier = Compute< false, false,
                        false, true,
                        false >(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 3:
         ier = Compute< false, false,
                        false, true,
                        true >(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 4:
         ier = Compute< false, false,
                        true, false,
                        false >(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 5:
         ier = Compute< false, false,
                        true, false,
                        true >(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 6:
         ier = Compute< false, false,
                        true, true,
                        false >(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 7:
         ier = Compute< false, false,
                        true, true,
                        true >(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 8:
         ier = Compute< false, true,
                        false, false,
                        false >(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 9:
         ier = Compute< false, true,
                        false, false,
                        true >(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 10:
         ier = Compute< false, true,
                        false, true,
                        false >(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 11:
         ier = Compute< false, true,
                        false, true,
                        true >(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 12:
         ier = Compute< false, true,
                        true, false,
                        false >(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 13:
         ier = Compute< false, true,
                        true, false,
                        true >(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 14:
         ier = Compute< false, true,
                        true, true,
                        false >(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 15:
         ier = Compute< false, true,
                        true, true,
                        true >(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 16:
         ier = Compute< true, false,
                        false, false,
                        false >(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 17:
         ier = Compute< true, false,
                        false, false,
                        true >(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 18:
         ier = Compute< true, false,
                        false, true,
                        false >(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 19:
         ier = Compute< true, false,
                        false, true,
                        true >(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 20:
         ier = Compute< true, false,
                        true, false,
                        false >(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 21:
         ier = Compute< true, false,
                        true, false,
                        true >(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 22:
         ier = Compute< true, false,
                        true, true,
                        false >(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 23:
         ier = Compute< true, false,
                        true, true,
                        true >(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 24:
         ier = Compute< true, true,
                        false, false,
                        false >(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 25:
         ier = Compute< true, true,
                        false, false,
                        true >(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 26:
         ier = Compute< true, true,
                        false, true,
                        false >(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 27:
         ier = Compute< true, true,
                        false, true,
                        true >(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 28:
         ier = Compute< true, true,
                        true, false,
                        false >(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 29:
         ier = Compute< true, true,
                        true, false,
                        true >(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 30:
         ier = Compute< true, true,
                        true, true,
                        false >(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 31:
         ier = Compute< true, true,
                        true, true,
                        true >(
                  modelCompute,
                  modelComputeArguments,
                  particleSpeciesCodes,
                  particleContributing,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      default:
         std::cout << "Unknown compute function index" << std::endl;
         ier = true;
         break;
   }
