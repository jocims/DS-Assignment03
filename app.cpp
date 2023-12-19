#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>

int inverseTable[] = {1, 9, 21, 15, 3, 19, 7, 23, 11, 5, 17, 25};
std::string ciphertext = "Dfanoryhmduy";
int *inverse = new int[4];
int **blocks;
int *blocksPerNode;
int *displacements;
int blockCount;
int *decodedBlocks;
int *ciphertextArray;
int numBlocks;
int totalDecodedBlocks = 0;

// printMatrix() will print the contents of the 2x2 matrix.
void printMatrix(int *array)
{
    std::cout << array[0] << " " << array[1] << std::endl;
    std::cout << array[2] << " " << array[3] << std::endl;
}

// printArray() will print the contents of an array of integers.
void printArray(int *arr, int size)
{
    for (int i = 0; i < size; i++)
    {
        std::cout << arr[i] << " ";
    }
    std::cout << std::endl;
}

// decodeBlock() will multiply a block of two integers by an inverse matrix, mod the result by 26, and return a decoded block of ints.
int *decodeBlock(int *block, int *inverse)
{
    int *decodedBlock = new int[2];

    decodedBlock[0] = ((inverse[0] * block[0] + inverse[1] * block[1])) % 26;
    decodedBlock[1] = ((inverse[2] * block[0] + inverse[3] * block[1])) % 26;

    return decodedBlock;
}

// calculateInvDeterminant() will calculate the modulo 26 inverse determinant of a 2x2 matrix.
int calculateInvDeterminant(int *matrix)
{
    // Calculate the determinant of the matrix
    int determinant = (matrix[0] * matrix[3] - matrix[1] * matrix[2] % 26);

    // Guarantee positive result
    if (determinant < 0)
    {
        determinant += 26;
    }

    int invDeterminant = -1;

    for (int i = 0; i < 12; i++)
    {
        if ((determinant * inverseTable[i]) % 26 == 1)
        {
            invDeterminant = inverseTable[i];
            break;
        }
    }
    return invDeterminant;
}

// calculateInverse() will calculate the modulo 26 inverse of a 2x2 matrix.
int *calculateInverse(int *matrix, int invDeterminant)
{

    int *inverse = new int[4];
    inverse[0] = (matrix[3] * invDeterminant) % 26;
    inverse[1] = (-matrix[1] * invDeterminant % 26) + 26; // Ensure positive result
    inverse[2] = (-matrix[2] * invDeterminant % 26) + 26; // Ensure positive result
    inverse[3] = (matrix[0] * invDeterminant) % 26;

    return inverse;
}

// Calculate the identity matrix
int *calculateIdentity(int *matrix, int *inverse)
{
    int *identity = new int[4];
    identity[0] = (matrix[0] * inverse[0] + matrix[1] * inverse[2]) % 26;
    identity[1] = (matrix[0] * inverse[1] + matrix[1] * inverse[3]) % 26;
    identity[2] = (matrix[2] * inverse[0] + matrix[3] * inverse[2]) % 26;
    identity[3] = (matrix[2] * inverse[1] + matrix[3] * inverse[3]) % 26;

    return identity;
}

int main(int argc, char **argv)
{
    // Initialize MPI
    MPI_Init(&argc, &argv);

    // Get the world rank and size
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Designate rank 0 as the coordinator and all other ranks as participants
    if (world_rank == 0)
    {
        // a)	Read in a matrix from the user, and use the helper method calculateInvDeterminant() to find its mod 26 inverse determinant.
        int matrix[4];
        std::cout << "Enter the 4 numbers of the matrix: ";
        std::cin >> matrix[0] >> matrix[1] >> matrix[2] >> matrix[3];
        int invDeterminant = calculateInvDeterminant(matrix);

        // print the inverse determinant
        std::cout << "The inverse determinant is: " << invDeterminant << std::endl;

        // use the helper method calculateInverse() to calculate the modulo 26 matrix inverse (A-1).
        inverse = calculateInverse(matrix, invDeterminant);

        // Calculate A x Inverse of A = Identity Matrix
        int *identity = calculateIdentity(matrix, inverse);

        // print the identity matrix
        std::cout << "The identity matrix is: " << std::endl;
        printMatrix(identity);

        // Check if the result is the identity matrix
        bool isIdentity = (identity[0] == 1 && identity[1] == 0 && identity[2] == 0 && identity[3] == 1);

        // Print the inverse matrix
        std::cout << "The inverse matrix is: " << std::endl;
        printMatrix(inverse);

        // The coordinator will then broadcast the inverse matrix to all other nodes.
        MPI_Bcast(inverse, 4, MPI_INT, 0, MPI_COMM_WORLD);

        // Convert the ciphertext to an array of integers using the provided mapping
        ciphertextArray = new int[ciphertext.length()];
        for (int i = 0; i < ciphertext.length(); i++)
        {
            ciphertextArray[i] = toupper(ciphertext[i]) - 'A'; // Direct mapping based on provided equivalents
        }

        printArray(ciphertextArray, ciphertext.length());

        // Split the ciphertext into blocks of size 2.
        numBlocks = ciphertext.length() / 2;
        blocks = new int *[numBlocks];
        for (int i = 0; i < numBlocks; i++)
        {
            // Allocate memory for each block (each block has 2 elements).
            blocks[i] = new int[2];

            // Assign values to each element in the block.
            blocks[i][0] = ciphertextArray[2 * i];
            blocks[i][1] = ciphertextArray[2 * i + 1];
        }

        int *decodedBlocks = new int[2 * numBlocks];

        // Calculate how many blocks node will receive.
        int blocksPerNode = numBlocks / world_size;
        int remainingBlocks = numBlocks % world_size;

        // Broadcast the value of remainingBlocks to all nodes
        MPI_Bcast(&remainingBlocks, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // Distribute the blocks to all nodes, including the coordinator.
        int currentBlockIndex = 0;
        int numBlocksToSend;

        for (int i = 0; i < world_size; i++)
        {
            // Calculate how many blocks to send to the current node.
            numBlocksToSend = (i < remainingBlocks) ? blocksPerNode + 1 : blocksPerNode;

            // Send the number of blocks to the current node.
            MPI_Send(&numBlocksToSend, 1, MPI_INT, i, 0, MPI_COMM_WORLD);

            // Send the actual blocks to the current node.
            for (int j = 0; j < numBlocksToSend; j++)
            {
                MPI_Send(&blocks[currentBlockIndex + j][0], 2, MPI_INT, i, 0, MPI_COMM_WORLD);
            }

            // Move to the next set of blocks.
            currentBlockIndex += numBlocksToSend;
        }

        // print only the blocks allocated for the coordinator
        for (int i = 0; i < blocksPerNode + 1; i++)
        {
            int *decodedBlock = decodeBlock(blocks[i], inverse);

            // Store the decoded blocks for node 0 in the final array.
            for (int j = 0; j < 2; j++)
            {
                decodedBlocks[totalDecodedBlocks + j] = decodedBlock[j];
            }
            totalDecodedBlocks += 2;
        }

        // Receive decoded blocks from participants.
        for (int i = 1; i < world_size; i++)
        {
            // Receive the number of decoded blocks.
            int numDecodedBlocks;
            MPI_Recv(&numDecodedBlocks, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Allocate memory for received decoded blocks.
            int *receivedDecodedBlocks = new int[2 * numDecodedBlocks];

            // Receive the decoded blocks.
            for (int j = 0; j < 2 * numDecodedBlocks; j++)
            {
                MPI_Recv(&receivedDecodedBlocks[j], 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            // Store the received decoded blocks in the final array.
            for (int j = 0; j < 2 * numDecodedBlocks; j++)
            {
                decodedBlocks[totalDecodedBlocks + j] = receivedDecodedBlocks[j];
            }
            totalDecodedBlocks += 2 * numDecodedBlocks;

            // Send acknowledgment to the participant that its decoded blocks are received.
            MPI_Send(&numDecodedBlocks, 1, MPI_INT, i, 0, MPI_COMM_WORLD);

            // Deallocate memory for received decoded blocks.
            delete[] receivedDecodedBlocks;
        }

        // Print or use the final array of decoded blocks for node 0.
        std::cout << "Node 0 final decoded blocks: ";
        printArray(decodedBlocks, totalDecodedBlocks);

        // Convert the decoded numbers to an array of letters using the provided mapping
        std::string decodedMessage;
        decodedMessage.reserve(totalDecodedBlocks / 2); // Assuming each block corresponds to a letter

        for (int i = 0; i < totalDecodedBlocks; i += 2)
        {
            char letter1 = 'A' + decodedBlocks[i];     // Convert the first integer to a letter
            char letter2 = 'A' + decodedBlocks[i + 1]; // Convert the second integer to a letter

            decodedMessage.push_back(letter1);
            decodedMessage.push_back(letter2);
        }

        // Print the decoded message
        std::cout << "The decoded message is: " << decodedMessage << std::endl;

        // Deallocate memory for the final array of decoded blocks.
        delete[] decodedBlocks;
    }
    else
    {

        // Receive the inverse matrix from the coordinator.
        MPI_Bcast(inverse, 4, MPI_INT, 0, MPI_COMM_WORLD);

        // Receive the number of blocks.
        int numBlocksToReceive;
        MPI_Recv(&numBlocksToReceive, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Allocate contiguous memory for two blocks.
        int *receivedBlocks = new int[2 * 2]; // Allocate for two blocks

        for (int i = 0; i < numBlocksToReceive; ++i)
        {
            MPI_Recv(&receivedBlocks[2 * i], 2, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // Decode the participant's blocks.
        for (int i = 0; i < numBlocksToReceive; i++)
        {
            // Calculate the start index of the current block in the receivedBlocks array.
            int blockStartIndex = 2 * i;

            // std::cout << "Node " << world_rank << " decoded blocks: ";
            int *decodedBlock = decodeBlock(&receivedBlocks[blockStartIndex], inverse);

            // Send the number of decoded blocks to the coordinator (node 0).
            MPI_Send(&numBlocksToReceive, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

            // Send the decoded block to the coordinator (node 0).
            for (int j = 0; j < 2; j++)
            {
                MPI_Send(&decodedBlock[j], 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                // print out what's being sent
                // std::cout << "Node " << world_rank << " MPI SEND: " << decodedBlock[j] << std::endl;
            }
        }

        // Deallocate memory.
        delete[] receivedBlocks;
    }

    // Finalize MPI
    MPI_Finalize();

    // Return status of 0 to the OS
    return 0;
}